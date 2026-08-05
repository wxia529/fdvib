// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package hessian assembles the force-constant matrix from central finite
// differences and writes the QE-format .dynG dynamical matrix, mirroring
// analysis.cpp.
package hessian

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/qeinput"
	"github.com/wxia529/fdvib/internal/qeoutput"
	"github.com/wxia529/fdvib/internal/settings"
	"github.com/wxia529/fdvib/internal/units"
)

// norm computes the length of a vector, like norm() in analysis.cpp.
func norm(v config.Vec3) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// WriteDynG writes the QE Gamma dynamical-matrix file with the exact
// iostream column layout of write_dynG().
func WriteDynG(p string, q *qeinput.QEInput, h []float64) error {
	alatAng := norm(q.Cell[0])
	alatBohr := alatAng / units.BohrToAng
	var b strings.Builder
	b.WriteString("Dynamical matrix file\nfdvib finite-difference local/Gamma matrix\n")
	b.WriteString(config.Right(fmt.Sprintf("%d", q.Ntyp), 3))
	b.WriteString(config.Right(fmt.Sprintf("%d", q.Nat), 5))
	b.WriteString(config.Right("0", 4))
	b.WriteString(config.Right(config.FormatFixed(alatBohr, 7), 12))
	for i := 0; i < 5; i++ {
		b.WriteString(config.Right(config.FormatFixed(0.0, 7), 12))
	}
	b.WriteString("\nBasis vectors\n")
	for _, v := range q.Cell {
		b.WriteString("  ")
		for k := 0; k < 3; k++ {
			b.WriteString(config.Right(config.FormatFixed(v[k]/alatAng, 9), 15))
		}
		b.WriteString("\n")
	}
	for i := 0; i < q.Ntyp; i++ {
		b.WriteString(config.Right(fmt.Sprintf("%d", i+1), 12))
		b.WriteString("  '")
		b.WriteString(config.Left(q.Species[i].Symbol, 7))
		b.WriteString("' ")
		b.WriteString(config.Right(config.FormatFixed(units.AmuRy*q.Species[i].Mass, 12), 22))
		b.WriteString("\n")
	}
	for i := 0; i < q.Nat; i++ {
		b.WriteString(config.Right(fmt.Sprintf("%d", i+1), 5))
		b.WriteString(config.Right(fmt.Sprintf("%d", q.Atoms[i].Type), 5))
		for k := 0; k < 3; k++ {
			b.WriteString(config.Right(config.FormatFixed(q.Atoms[i].R[k]/alatAng, 10), 18))
		}
		b.WriteString("\n")
	}
	b.WriteString("\n     Dynamical  Matrix in cartesian axes\n\n     q = (    0.000000000   0.000000000   0.000000000 ) \n\n")
	n3 := 3 * q.Nat
	for a := 0; a < q.Nat; a++ {
		for b2 := 0; b2 < q.Nat; b2++ {
			b.WriteString(config.Right(fmt.Sprintf("%d", a+1), 5))
			b.WriteString(config.Right(fmt.Sprintf("%d", b2+1), 5))
			b.WriteString("\n")
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					b.WriteString(config.Right(config.FormatFixed(h[(3*a+i)*n3+3*b2+j], 10), 14))
					b.WriteString("   ")
					b.WriteString(config.Right(config.FormatFixed(0.0, 10), 12))
					b.WriteString("  ")
				}
				b.WriteString("\n")
			}
		}
	}
	return config.WriteText(p, b.String())
}

// completedForces locates the forces.dat of a completed displacement job via
// its state marker, like completed_forces().
func completedForces(s *settings.Settings, id string) (string, error) {
	marker := filepath.Join(s.Workdir, "state", id+".complete")
	text, err := config.ReadText(marker)
	if err != nil {
		return "", err
	}
	fields := strings.Fields(text)
	if len(fields) != 2 || !strings.HasPrefix(fields[0], id+"_") ||
		filepath.Base(fields[0]) != fields[0] {
		return "", fmt.Errorf("Invalid displacement completion snapshot: %s", marker)
	}
	return filepath.Join(s.Workdir, "calculations", fields[0], "forces.dat"), nil
}

// Analyze assembles the Hessian from completed displacement jobs and writes
// the dynG/dynmat.in/metadata.dat results, like analyze().
func Analyze(s *settings.Settings) error {
	scfPath := filepath.Join(s.Workdir, "scf.in.reference")
	q, err := qeinput.ParseQeInput(scfPath)
	if err != nil {
		return err
	}
	selected := s.Selected
	if s.SystemType == "gas" {
		selected = make([]int, q.Nat)
		for i := range selected {
			selected[i] = i + 1
		}
	}
	n3 := 3 * q.Nat
	h := make([]float64, n3*n3)
	deltaBohr := s.Displacement / units.BohrToAng
	for _, atom1 := range selected {
		for axis := 0; axis < 3; axis++ {
			pPath, err := completedForces(s, settings.JobName(atom1, axis, 1))
			if err != nil {
				return err
			}
			mPath, err := completedForces(s, settings.JobName(atom1, axis, -1))
			if err != nil {
				return err
			}
			p, err := qeoutput.ReadForces(pPath, q.Nat)
			if err != nil {
				return err
			}
			m, err := qeoutput.ReadForces(mPath, q.Nat)
			if err != nil {
				return err
			}
			col := 3*(atom1-1) + axis
			for _, atomj1 := range selected {
				for beta := 0; beta < 3; beta++ {
					row := 3*(atomj1-1) + beta
					h[row*n3+col] = -(p[atomj1-1][beta] - m[atomj1-1][beta]) / (2 * deltaBohr)
				}
			}
		}
	}
	asym := 0.0
	for i := 0; i < n3; i++ {
		for j := 0; j < i; j++ {
			if math.Abs(h[i*n3+j]-h[j*n3+i]) > asym {
				asym = math.Abs(h[i*n3+j] - h[j*n3+i])
			}
			x := 0.5 * (h[i*n3+j] + h[j*n3+i])
			h[i*n3+j] = x
			h[j*n3+i] = x
		}
	}
	resultsDir := filepath.Join(s.Workdir, "results")
	if err := os.MkdirAll(resultsDir, 0o755); err != nil {
		return fmt.Errorf("Cannot write %s", resultsDir)
	}
	dyn := filepath.Join(resultsDir, s.OutputPrefix+".dynG")
	for _, p := range []string{dyn, filepath.Join(resultsDir, "dynmat.in"),
		filepath.Join(resultsDir, "metadata.dat")} {
		if _, err := os.Stat(p); err == nil {
			return fmt.Errorf("Refuse to overwrite %s", p)
		}
	}
	metadataText, err := config.ReadText(filepath.Join(s.Workdir, "metadata.dat"))
	if err != nil {
		return err
	}
	// fs::copy_file preserves the source permission bits.
	if fi, statErr := os.Stat(filepath.Join(s.Workdir, "metadata.dat")); statErr == nil {
		if err := os.WriteFile(filepath.Join(resultsDir, "metadata.dat"), []byte(metadataText), fi.Mode().Perm()); err != nil {
			return fmt.Errorf("Cannot write %s", filepath.Join(resultsDir, "metadata.dat"))
		}
	} else {
		return statErr
	}
	if err := WriteDynG(dyn, q, h); err != nil {
		return err
	}
	remove := ".true."
	if s.SystemType == "gas" {
		remove = ".false."
	}
	di := fmt.Sprintf("&INPUT\n  fildyn='%s',\n  filout='%s.freq.out',\n  asr='no',\n  remove_interaction_blocks=%s,\n/\n",
		filepath.Base(dyn), s.OutputPrefix, remove)
	if err := config.WriteText(filepath.Join(resultsDir, "dynmat.in"), di); err != nil {
		return err
	}
	fmt.Printf("Hessian max asymmetry: %s Ry/Bohr^2\n", config.FormatSci(asym, 6))
	fmt.Printf("Wrote %s and %s\n", config.DisplayPath(dyn),
		config.DisplayPath(filepath.Join(resultsDir, "dynmat.in")))
	return nil
}
