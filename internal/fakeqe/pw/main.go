// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Command pw is a fake pw.x: it evaluates forces from a fixed quadratic
// potential and prints a QE-style output.
package main

import (
	"fmt"
	"math"
	"os"

	"github.com/wxia529/fdvib/internal/fakeqe"
)

func main() {
	// Arguments arrive as "fakeqe-pw -inp scf.in".
	input := ""
	for i := 1; i < len(os.Args)-1; i++ {
		if os.Args[i] == "-inp" {
			input = os.Args[i+1]
			break
		}
	}
	if input == "" {
		fmt.Fprintln(os.Stderr, "fakeqe-pw: missing -inp argument")
		os.Exit(127)
	}
	symbols, types, coords, prefix, err := fakeqe.ParseInput(input)
	if err != nil {
		fmt.Fprintln(os.Stderr, "fakeqe-pw:", err)
		os.Exit(127)
	}
	nat := len(coords)
	ref, err := fakeqe.LoadReference()
	if err != nil {
		// First run (reference SCF): record the coordinates.
		if err := fakeqe.SaveReference(coords); err != nil {
			fmt.Fprintln(os.Stderr, "fakeqe-pw:", err)
			os.Exit(127)
		}
		ref = coords
	}
	if len(ref) != nat {
		fmt.Fprintln(os.Stderr, "fakeqe-pw: reference coordinate count mismatch")
		os.Exit(1)
	}
	if err := fakeqe.WriteDensity(".", prefix); err != nil {
		fmt.Fprintln(os.Stderr, "fakeqe-pw:", err)
		os.Exit(1)
	}

	h := fakeqe.Hessian(nat)
	n3 := 3 * nat
	forces := make([][3]float64, nat)
	for i := 0; i < nat; i++ {
		for j := 0; j < nat; j++ {
			for alpha := 0; alpha < 3; alpha++ {
				for beta := 0; beta < 3; beta++ {
					// Coordinates are in Angstrom; H is in Ry/Bohr^2.
					delta := (coords[j][beta] - ref[j][beta]) / fakeqe.BohrToAng
					forces[i][alpha] -= h[(3*i+alpha)*n3+3*j+beta] * delta
				}
			}
		}
	}

	fmt.Println("     Program PWSCF v.7.2")
	fmt.Println()
	fmt.Println("     Self-consistency achieved")
	fmt.Println()
	fmt.Println("     Forces acting on atoms (cartesian axes, Ry/au):")
	fmt.Println()
	total := 0.0
	for i := 0; i < nat; i++ {
		t := 0
		if i < len(types) {
			t = types[i]
		}
		_ = symbols
		fmt.Printf("     atom %4d type %2d   force = %16.10f %16.10f %16.10f\n",
			i+1, t, forces[i][0], forces[i][1], forces[i][2])
		total += forces[i][0]*forces[i][0] + forces[i][1]*forces[i][1] + forces[i][2]*forces[i][2]
	}
	fmt.Println()
	fmt.Printf("     Total force = %14.8f     Total stress = %14.8f\n\n", math.Sqrt(total), 0.0)
	fmt.Println("!    total energy              =     -17.12345678901234 Ry")
	fmt.Println()
	fmt.Println("     JOB DONE.")
	os.Exit(0)
}
