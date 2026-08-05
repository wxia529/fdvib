// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Command dynmat is a fake dynmat.x: it mass-weights and diagonalizes the
// .dynG matrix and writes QE-style frequency output.
package main

import (
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mat"

	"github.com/wxia529/fdvib/internal/fakeqe"
)

type dynGeometry struct {
	nat    int
	masses []float64 // amu
	matrix []float64 // 3N x 3N cartesian, Ry/Bohr^2
}

// parseDynG reads a .dynG file written by FDVIB (or QE).
func parseDynG(path string) (*dynGeometry, error) {
	text, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(string(text), "\n")
	if len(lines) < 4 {
		return nil, fmt.Errorf("short dynG")
	}
	header := strings.Fields(lines[2])
	if len(header) < 4 {
		return nil, fmt.Errorf("bad dynG header")
	}
	ntyp, err1 := strconv.Atoi(header[0])
	nat, err2 := strconv.Atoi(header[1])
	alat, err3 := strconv.ParseFloat(header[3], 64)
	if err1 != nil || err2 != nil || err3 != nil {
		return nil, fmt.Errorf("bad dynG header numbers")
	}
	_ = ntyp
	typeMasses := make([]float64, 0)
	for i := 0; i < ntyp; i++ {
		fields := strings.Fields(lines[7+i])
		if len(fields) == 0 {
			return nil, fmt.Errorf("bad species row")
		}
		massRy, err := strconv.ParseFloat(fields[len(fields)-1], 64)
		if err != nil {
			return nil, err
		}
		typeMasses = append(typeMasses, massRy/911.4442431390707)
	}
	_ = alat
	g := &dynGeometry{nat: nat}
	g.masses = make([]float64, nat)
	for i := 0; i < nat; i++ {
		fields := strings.Fields(lines[7+ntyp+i])
		if len(fields) < 5 {
			return nil, fmt.Errorf("bad atom row")
		}
		ti, _ := strconv.Atoi(fields[1])
		g.masses[i] = typeMasses[ti-1]
	}
	// Dynamical matrix block: "     Dynamical  Matrix in cartesian axes".
	blockStart := -1
	for i, line := range lines {
		if strings.Contains(line, "Dynamical  Matrix in cartesian axes") {
			blockStart = i
			break
		}
	}
	if blockStart < 0 {
		return nil, fmt.Errorf("no matrix block in dynG")
	}
	n3 := 3 * nat
	g.matrix = make([]float64, n3*n3)
	row := 0
	col := 0
	started := false
	for i := blockStart; i < len(lines); i++ {
		fields := strings.Fields(lines[i])
		if len(fields) == 2 && started {
			// Pair header "a b" starts a new 3x3 block.
			a, _ := strconv.Atoi(fields[0])
			b, _ := strconv.Atoi(fields[1])
			row = 3 * (a - 1)
			col = 3 * (b - 1)
			continue
		}
		if !started {
			// First pair header begins the matrix.
			if len(fields) == 2 {
				started = true
				a, _ := strconv.Atoi(fields[0])
				b, _ := strconv.Atoi(fields[1])
				row = 3 * (a - 1)
				col = 3 * (b - 1)
			}
			continue
		}
		if len(fields) >= 6 {
			// Three (re im) pairs for one cartesian row.
			for k := 0; k < 3; k++ {
				re, err := strconv.ParseFloat(fields[2*k], 64)
				if err != nil {
					return nil, err
				}
				g.matrix[row*n3+col+k] = re
			}
			row++
		}
	}
	return g, nil
}

func main() {
	input := ""
	for i := 1; i < len(os.Args)-1; i++ {
		if os.Args[i] == "-inp" {
			input = os.Args[i+1]
			break
		}
	}
	if input == "" {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat: missing -inp argument")
		os.Exit(127)
	}
	// Read dynmat.in to find fildyn/filout.
	text, err := os.ReadFile(input)
	if err != nil {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat:", err)
		os.Exit(127)
	}
	fildyn, filout := "", ""
	for _, line := range strings.Split(string(text), "\n") {
		line = strings.TrimSpace(strings.TrimSuffix(line, ","))
		fields := strings.SplitN(line, "=", 2)
		if len(fields) != 2 {
			continue
		}
		key := strings.ToLower(strings.TrimSpace(fields[0]))
		val := strings.Trim(strings.TrimSpace(fields[1]), "'\"")
		if key == "fildyn" {
			fildyn = val
		} else if key == "filout" {
			filout = val
		}
	}
	if fildyn == "" || filout == "" {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat: missing fildyn/filout in dynmat.in")
		os.Exit(1)
	}
	g, err := parseDynG(fildyn)
	if err != nil {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat:", err)
		os.Exit(1)
	}

	// Mass-weighted dynamical matrix: D = M^{-1/2} C M^{-1/2}.
	n3 := 3 * g.nat
	d := mat.NewSymDense(n3, nil)
	for i := 0; i < n3; i++ {
		for j := 0; j <= i; j++ {
			val := g.matrix[i*n3+j] / math.Sqrt(g.masses[i/3]*g.masses[j/3])
			d.SetSym(i, j, val)
		}
	}
	var es mat.EigenSym
	ok := es.Factorize(d, true)
	if !ok {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat: eigen factorization failed")
		os.Exit(1)
	}
	values := es.Values(nil)
	var vectors mat.Dense
	es.VectorsTo(&vectors)

	// Frequency conversion: omega^2 = lambda * (Ry/Bohr^2 -> SI) / (mass in
	// Ry units -> kg); nu = omega / 2pi.
	ryPerBohr2 := fakeqe.RyToJ / (fakeqe.BohrToM * fakeqe.BohrToM)
	var out strings.Builder
	out.WriteString("     Diagonalizing the dynamical matrix\n\n")
	for m := 0; m < n3; m++ {
		omega2 := values[m] * ryPerBohr2 / fakeqe.MeKg
		nuHz := math.Sqrt(math.Abs(omega2)) / (2 * math.Pi)
		sign := 1.0
		if values[m] < 0 {
			sign = -1.0
		}
		thz := sign * nuHz / 1e12
		cm1 := sign * nuHz / fakeqe.CCmPerS
		fmt.Fprintf(&out, "   freq (%4d) = %10.6f [THz] = %10.6f [cm-1]\n", m+1, thz, cm1)
		for a := 0; a < g.nat; a++ {
			x := vectors.At(3*a, m)
			y := vectors.At(3*a+1, m)
			z := vectors.At(3*a+2, m)
			fmt.Fprintf(&out, "( %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f )\n",
				x, 0.0, y, 0.0, z, 0.0)
		}
	}
	out.WriteString("\n     JOB DONE.\n")
	// dynmat.out is captured by FDVIB; freq.out is written next to it.
	if err := os.WriteFile("dynmat.out", []byte(out.String()), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat:", err)
		os.Exit(1)
	}
	if err := os.WriteFile(filout, []byte(out.String()), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, "fakeqe-dynmat:", err)
		os.Exit(1)
	}
	os.Exit(0)
}
