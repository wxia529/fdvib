// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package modes implements rigid-body classification and normal-mode
// selection for SHM/Molden export, mirroring mode_selection.cpp. The
// inertia eigenproblems use the same stable analytic 3x3 closed form as the
// C++ implementation so classification results match exactly.
package modes

import (
	"fmt"
	"math"
	"sort"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/qevibration"
	"github.com/wxia529/fdvib/internal/results"
	"github.com/wxia529/fdvib/internal/units"
)

// ModeSelection is the outcome of a mode selection, like struct ModeSelection.
type ModeSelection struct {
	Indices        []int
	Classification string
	LargestRemoved float64
}

// inertiaEigenvalues computes the sorted principal moments of inertia of a
// geometry (Bohr/amu) about its center of mass, like inertia_eigenvalues().
func inertiaEigenvalues(g *qevibration.DynGeometry, context string) ([3]float64, error) {
	totalMass := 0.0
	for _, m := range g.Masses {
		totalMass += m
	}
	if !(totalMass > 0.0) {
		return [3]float64{}, fmt.Errorf("invalid total mass for %s", context)
	}
	var center config.Vec3
	for i := range g.Masses {
		for k := 0; k < 3; k++ {
			center[k] += g.Masses[i] * g.RBohr[i][k] / totalMass
		}
	}
	var diagonal, off config.Vec3
	for i := range g.Masses {
		x := g.RBohr[i][0] - center[0]
		y := g.RBohr[i][1] - center[1]
		z := g.RBohr[i][2] - center[2]
		mass := g.Masses[i]
		diagonal[0] += mass * (y*y + z*z)
		diagonal[1] += mass * (x*x + z*z)
		diagonal[2] += mass * (x*x + y*y)
		off[0] -= mass * x * y
		off[1] -= mass * x * z
		off[2] -= mass * y * z
	}
	return symmetricEigenvalues(diagonal, off), nil
}

// symmetricEigenvalues is the stable analytic eigenvalue solution for a real
// symmetric 3x3 matrix given as diagonal and off-diagonal components,
// matching symmetric_eigenvalues() in thermo.cpp and inertia_eigenvalues()
// in mode_selection.cpp.
func symmetricEigenvalues(diag, off config.Vec3) [3]float64 {
	p1 := off[0]*off[0] + off[1]*off[1] + off[2]*off[2]
	if p1 == 0 {
		a := [3]float64{diag[0], diag[1], diag[2]}
		sort.Float64s(a[:])
		return a
	}
	q := (diag[0] + diag[1] + diag[2]) / 3
	p2 := (diag[0]-q)*(diag[0]-q) + (diag[1]-q)*(diag[1]-q) +
		(diag[2]-q)*(diag[2]-q) + 2.0*p1
	p := math.Sqrt(p2 / 6.0)
	a00 := (diag[0] - q) / p
	a11 := (diag[1] - q) / p
	a22 := (diag[2] - q) / p
	a01 := off[0] / p
	a02 := off[1] / p
	a12 := off[2] / p
	det := a00*a11*a22 + 2*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01
	phi := math.Acos(clamp(det/2.0, -1.0, 1.0)) / 3
	e3 := q + 2*p*math.Cos(phi)
	e1 := q + 2*p*math.Cos(phi+2*units.PI/3)
	e2 := 3*q - e1 - e3
	a := [3]float64{e1, e2, e3}
	sort.Float64s(a[:])
	return a
}

func clamp(x, lo, hi float64) float64 {
	if x < lo {
		return lo
	}
	if x > hi {
		return hi
	}
	return x
}

// retainLargestModes keeps the modes with the largest absolute frequency,
// like retain_largest_modes().
func retainLargestModes(modes []qevibration.Mode, keep int, classification, context string) (*ModeSelection, error) {
	if keep > len(modes) {
		return nil, fmt.Errorf("requested %s vibration count exceeds available modes", context)
	}
	order := make([]int, len(modes))
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(a, b int) bool {
		return math.Abs(modes[order[a]].Freq) < math.Abs(modes[order[b]].Freq)
	})
	removed := make([]bool, len(modes))
	selected := &ModeSelection{Classification: classification}
	for i := 0; i < len(modes)-keep; i++ {
		removed[order[i]] = true
		if math.Abs(modes[order[i]].Freq) > selected.LargestRemoved {
			selected.LargestRemoved = math.Abs(modes[order[i]].Freq)
		}
	}
	for i := range modes {
		if !removed[i] {
			selected.Indices = append(selected.Indices, i)
		}
	}
	return selected, nil
}

// selectNonzeroModes keeps modes with nonzero frequency, like
// select_nonzero_modes().
func selectNonzeroModes(modes []qevibration.Mode, classification string) *ModeSelection {
	selected := &ModeSelection{Classification: classification}
	for i := range modes {
		if modes[i].Freq != 0.0 {
			selected.Indices = append(selected.Indices, i)
		}
	}
	return selected
}

// metadataSelectedAtoms parses the selected atom list from metadata.dat,
// like metadata_selected_atoms().
func metadataSelectedAtoms(metadata *results.ResultMetadata, nat int) ([]int, error) {
	var selected []int
	if metadata.SelectedAtoms == "all" {
		selected = make([]int, nat)
		for i := range selected {
			selected[i] = i + 1
		}
	} else {
		var err error
		selected, err = config.IntegerList(metadata.SelectedAtoms)
		if err != nil {
			return nil, err
		}
	}
	if len(selected) == 0 {
		return nil, fmt.Errorf("local SHM export requires selected_atoms in metadata.dat")
	}
	unique := make(map[int]bool)
	for _, atom := range selected {
		if atom < 1 || atom > nat || unique[atom] {
			return nil, fmt.Errorf("invalid or duplicate selected atom in metadata.dat")
		}
		unique[atom] = true
	}
	return selected, nil
}

// SelectGasInternalModes selects the internal (rigid-body-free) modes for a
// gas-phase molecule, like select_gas_internal_modes().
func SelectGasInternalModes(geometry *qevibration.DynGeometry, modes []qevibration.Mode,
	context string) (*ModeSelection, error) {
	nat := len(geometry.Masses)
	if nat == 1 {
		return &ModeSelection{Classification: "atom"}, nil
	}
	inertia, err := inertiaEigenvalues(geometry, context)
	if err != nil {
		return nil, err
	}
	linear := false
	for _, value := range inertia {
		if value < 0.001 {
			linear = true
			break
		}
	}
	classification := "nonlinear"
	keep := 3*nat - 6
	if linear {
		classification = "linear"
		keep = 3*nat - 5
	}
	return retainLargestModes(modes, keep, classification, context)
}

// SelectShmModes selects modes for Shermo export, like select_shm_modes().
func SelectShmModes(metadata *results.ResultMetadata, geometry *qevibration.DynGeometry,
	modes []qevibration.Mode) (*ModeSelection, error) {
	switch metadata.ModeSelection {
	case "gas":
		return SelectGasInternalModes(geometry, modes, "SHM export")
	case "local":
		selected, err := metadataSelectedAtoms(metadata, len(geometry.Masses))
		if err != nil {
			return nil, err
		}
		return retainLargestModes(modes, 3*len(selected), "local", "SHM export")
	case "all":
		return selectNonzeroModes(modes, "all"), nil
	default:
		return nil, fmt.Errorf("mode_selection must be all, gas, or local in metadata.dat")
	}
}

// SelectMoldenModes selects modes for Molden export, like
// select_molden_modes().
func SelectMoldenModes(metadata *results.ResultMetadata, geometry *qevibration.DynGeometry,
	modes []qevibration.Mode) (*ModeSelection, error) {
	switch metadata.ModeSelection {
	case "gas":
		return SelectGasInternalModes(geometry, modes, "Molden export")
	case "local", "all":
		return selectNonzeroModes(modes, metadata.ModeSelection), nil
	default:
		return nil, fmt.Errorf("mode_selection must be all, gas, or local in metadata.dat")
	}
}

// SymmetricEigenvalues is the exported closed-form 3x3 symmetric eigen
// solver (used by the thermochemistry package), like symmetric_eigenvalues()
// in thermo.cpp.
func SymmetricEigenvalues(diag, off config.Vec3) [3]float64 {
	return symmetricEigenvalues(diag, off)
}
