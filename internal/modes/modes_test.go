// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package modes

import (
	"testing"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/qevibration"
	"github.com/wxia529/fdvib/internal/results"
)

func freqMode(freqs ...float64) []qevibration.Mode {
	modes := make([]qevibration.Mode, len(freqs))
	for i, f := range freqs {
		modes[i].Freq = f
		modes[i].Displacement = []config.Vec3{{0, 0, 0}}
	}
	return modes
}

// TestSymmetricEigenvalues checks the analytic 3x3 eigen solver against
// known values (diagonal matrix -> sorted diagonal).
func TestSymmetricEigenvalues(t *testing.T) {
	got := SymmetricEigenvalues([3]float64{3, 1, 2}, [3]float64{0, 0, 0})
	if got[0] != 1 || got[1] != 2 || got[2] != 3 {
		t.Errorf("diagonal eigen = %v", got)
	}
	// Off-diagonal case with a known spectrum.
	// Matrix [[2,1,0],[1,2,1],[0,1,2]] has eigenvalues 2-sqrt(2), 2, 2+sqrt(2).
	got = SymmetricEigenvalues([3]float64{2, 2, 2}, [3]float64{-1, 0, -1})
	want := [3]float64{2 - 1.4142135623730951, 2, 2 + 1.4142135623730951}
	for i := range got {
		if abs(got[i]-want[i]) > 1e-12 {
			t.Errorf("eigen[%d] = %v, want %v", i, got, want)
		}
	}
}

func TestSelectGasInternalModes(t *testing.T) {
	// Linear CO2: 3 atoms -> 3N-5 = 4 internal modes.
	geo := &qevibration.DynGeometry{
		Masses: []float64{12.0, 16.0, 16.0},
		RBohr:  []config.Vec3{{0, 0, 0}, {2, 0, 0}, {4, 0, 0}},
	}
	modes := freqMode(0, 0, 0, 0, 0, 100, 300, 700, 1400)
	sel, err := SelectGasInternalModes(geo, modes, "test")
	if err != nil {
		t.Fatal(err)
	}
	if sel.Classification != "linear" || len(sel.Indices) != 4 {
		t.Errorf("selection: %+v", sel)
	}
	// Largest removed |freq| is 0 (all zero modes removed first).
	if sel.LargestRemoved != 0 {
		t.Errorf("largest removed = %v", sel.LargestRemoved)
	}
	if sel.Indices[0] != 5 || sel.Indices[3] != 8 {
		t.Errorf("indices = %v", sel.Indices)
	}

	// Nonlinear H2O-like: 3N-6 = 3.
	geo = &qevibration.DynGeometry{
		Masses: []float64{16.0, 1.0, 1.0},
		RBohr:  []config.Vec3{{0, 0, 0}, {2, 0, 0}, {0, 2, 0}},
	}
	sel, err = SelectGasInternalModes(geo, freqMode(0, 0, 0, 0, 0, 0, 100, 300, 700), "test")
	if err != nil {
		t.Fatal(err)
	}
	if sel.Classification != "nonlinear" || len(sel.Indices) != 3 {
		t.Errorf("selection: %+v", sel)
	}

	// Monatomic: classification atom, no modes.
	geo = &qevibration.DynGeometry{Masses: []float64{40.0}}
	sel, err = SelectGasInternalModes(geo, freqMode(0, 0, 0), "test")
	if err != nil {
		t.Fatal(err)
	}
	if sel.Classification != "atom" || len(sel.Indices) != 0 {
		t.Errorf("selection: %+v", sel)
	}
}

func TestSelectShmModes(t *testing.T) {
	geo := &qevibration.DynGeometry{
		Masses: []float64{12.0, 16.0, 16.0},
		RBohr:  []config.Vec3{{0, 0, 0}, {2, 0, 0}, {4, 0, 0}},
	}
	modes := freqMode(0, 0, 0, 0, 0, 100, 300, 700, 1400)

	// local with selected_atoms 1,2: keep 6 modes.
	md := &results.ResultMetadata{Program: "qe", Multiplicity: 1, ModeSelection: "local", SelectedAtoms: "1,2"}
	sel, err := SelectShmModes(md, geo, modes)
	if err != nil {
		t.Fatal(err)
	}
	if sel.Classification != "local" || len(sel.Indices) != 6 {
		t.Errorf("selection: %+v", sel)
	}
	if sel.LargestRemoved != 0 {
		t.Errorf("largest removed = %v", sel.LargestRemoved)
	}

	// all: nonzero only.
	md = &results.ResultMetadata{Program: "qe", Multiplicity: 1, ModeSelection: "all", SelectedAtoms: "all"}
	sel, err = SelectShmModes(md, geo, modes)
	if err != nil {
		t.Fatal(err)
	}
	if sel.Classification != "all" || len(sel.Indices) != 4 {
		t.Errorf("selection: %+v", sel)
	}

	// local without selected_atoms fails.
	md = &results.ResultMetadata{Program: "qe", Multiplicity: 1, ModeSelection: "local", SelectedAtoms: ""}
	if _, err := SelectShmModes(md, geo, modes); err == nil {
		t.Error("local without selected_atoms should fail")
	}
}

func TestSelectMoldenModes(t *testing.T) {
	geo := &qevibration.DynGeometry{Masses: []float64{16.0, 1.0, 1.0}}
	modes := freqMode(0, 0, 0, 0, 0, 0, 100, 300, 700)
	md := &results.ResultMetadata{Program: "qe", Multiplicity: 1, ModeSelection: "all", SelectedAtoms: "all"}
	sel, err := SelectMoldenModes(md, geo, modes)
	if err != nil {
		t.Fatal(err)
	}
	if len(sel.Indices) != 3 {
		t.Errorf("selection: %+v", sel)
	}
	md = &results.ResultMetadata{Program: "qe", Multiplicity: 1, ModeSelection: "bogus", SelectedAtoms: "all"}
	if _, err := SelectMoldenModes(md, geo, modes); err == nil {
		t.Error("bogus mode_selection should fail")
	}
}

func TestRetainLargestFailsWhenTooFew(t *testing.T) {
	geo := &qevibration.DynGeometry{Masses: []float64{16.0, 1.0, 1.0}}
	md := &results.ResultMetadata{Program: "qe", Multiplicity: 1, ModeSelection: "local", SelectedAtoms: "1,2,3,4"}
	if _, err := SelectShmModes(md, geo, freqMode(0, 0, 0)); err == nil {
		t.Error("keep > available should fail")
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
