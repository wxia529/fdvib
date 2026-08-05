// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package fakeqe is a test double for Quantum ESPRESSO. The pw subcommand
// evaluates forces from a fixed quadratic potential F = -H*(r - r_ref) so
// the finite-difference Hessian recovered by FDVIB is exactly H; the dynmat
// subcommand mass-weights and diagonalizes the .dynG matrix with gonum to
// emit a QE-style frequency output. Both are driven through the FDVIB
// pw_command/dynmat_command settings.
package fakeqe

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
)

// Physical constants used by the dynmat frequency conversion (CODATA).
const (
	RyToJ     = 4.3597447222071e-18 // J per Ry
	BohrToM   = 5.29177210903e-11   // m per Bohr
	BohrToAng = 0.529177210544      // angstrom per bohr
	MeKg      = 9.1093837015e-31    // electron mass in kg (1 Ry mass unit)
	CCmPerS   = 2.99792458e10       // speed of light in cm/s
)

// StateEnv is the environment variable naming the fake-QE state directory.
const StateEnv = "FDVIB_FAKE_STATE"

// HEnvA/HEnvB/HEnvC tune the fake quadratic potential:
// H[(3i+a),(3j+b)] = a*delta(i,j)*delta(a,b) + b*delta(i,j) +
// c*delta(|i-j|,1), where delta(x,y) is the Kronecker delta (Ry/Bohr^2).
const (
	HEnvA = "FDVIB_FAKE_A"
	HEnvB = "FDVIB_FAKE_B"
	HEnvC = "FDVIB_FAKE_C"
)

var (
	speciesCardRe   = regexp.MustCompile(`(?i)^\s*ATOMIC_SPECIES\b`)
	positionsCardRe = regexp.MustCompile(`(?i)^\s*ATOMIC_POSITIONS\b`)
	prefixRe        = regexp.MustCompile(`(?i)\bprefix\s*=\s*['"]([^'"]+)['"]`)
)

func envFloat(name string, def float64) float64 {
	if v := os.Getenv(name); v != "" {
		if x, err := strconv.ParseFloat(v, 64); err == nil {
			return x
		}
	}
	return def
}

// ParseInput extracts species, coordinates, and prefix from a pw.x-style
// input file.
func ParseInput(path string) (symbols []string, types []int, coords [][3]float64, prefix string, err error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, nil, "", fmt.Errorf("cannot read input: %w", err)
	}
	defer f.Close()

	var lines []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	prefix = "pwscf"
	if m := prefixRe.FindStringSubmatch(strings.Join(lines, "\n")); m != nil {
		prefix = m[1]
	}
	var species []string
	posStart := -1
	for i, line := range lines {
		trimmed := strings.TrimSpace(line)
		if speciesCardRe.MatchString(trimmed) {
			species = nil
			for j := i + 1; j < len(lines); j++ {
				fields := strings.Fields(lines[j])
				if len(fields) == 0 || strings.HasPrefix(fields[0], "&") {
					break
				}
				if len(fields) >= 3 {
					species = append(species, fields[0])
				} else {
					break
				}
			}
		}
		if positionsCardRe.MatchString(trimmed) {
			posStart = i + 1
			break
		}
	}
	if posStart < 0 {
		return nil, nil, nil, "", fmt.Errorf("no ATOMIC_POSITIONS card")
	}
	for j := posStart; j < len(lines); j++ {
		fields := strings.Fields(lines[j])
		if len(fields) < 4 {
			break
		}
		symbols = append(symbols, fields[0])
		var c [3]float64
		for k := 0; k < 3; k++ {
			x, err := strconv.ParseFloat(fields[1+k], 64)
			if err != nil {
				return nil, nil, nil, "", fmt.Errorf("bad coordinate: %w", err)
			}
			c[k] = x
		}
		coords = append(coords, c)
	}
	for _, sym := range symbols {
		t := 0
		for i, s := range species {
			if s == sym {
				t = i + 1
				break
			}
		}
		types = append(types, t)
	}
	return symbols, types, coords, prefix, nil
}

// Hessian builds the fake force-constant matrix (Ry/Bohr^2) for nat atoms.
func Hessian(nat int) []float64 {
	a := envFloat(HEnvA, 0.4)
	b := envFloat(HEnvB, 0.05)
	c := envFloat(HEnvC, 0.02)
	n3 := 3 * nat
	h := make([]float64, n3*n3)
	for i := 0; i < nat; i++ {
		for j := 0; j < nat; j++ {
			for alpha := 0; alpha < 3; alpha++ {
				for beta := 0; beta < 3; beta++ {
					row := 3*i + alpha
					col := 3*j + beta
					var v float64
					if i == j {
						v = b
						if alpha == beta {
							v += a
						}
					} else if abs(i-j) == 1 {
						v = c
					}
					h[row*n3+col] = v
				}
			}
		}
	}
	return h
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// RefPath returns the reference-coordinate state file.
func RefPath() string {
	return filepath.Join(os.Getenv(StateEnv), "ref.xyz")
}

// LoadReference reads the saved reference coordinates.
func LoadReference() ([][3]float64, error) {
	f, err := os.Open(RefPath())
	if err != nil {
		return nil, fmt.Errorf("no fake reference state (run the reference SCF first): %w", err)
	}
	defer f.Close()
	var coords [][3]float64
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) != 3 {
			continue
		}
		var c [3]float64
		for k := 0; k < 3; k++ {
			x, err := strconv.ParseFloat(fields[k], 64)
			if err != nil {
				return nil, err
			}
			c[k] = x
		}
		coords = append(coords, c)
	}
	return coords, nil
}

// SaveReference stores the reference coordinates.
func SaveReference(coords [][3]float64) error {
	if dir := os.Getenv(StateEnv); dir != "" {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			return err
		}
	}
	var b strings.Builder
	for _, c := range coords {
		fmt.Fprintf(&b, "%.12f %.12f %.12f\n", c[0], c[1], c[2])
	}
	return os.WriteFile(RefPath(), []byte(b.String()), 0o644)
}

// WriteDensity writes a fake charge-density file (and optionally paw.txt
// when FDVIB_FAKE_PAW is set) for an attempt.
func WriteDensity(cwd, prefix string) error {
	p := filepath.Join(cwd, "out", prefix+".save", "charge-density.dat")
	if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
		return err
	}
	if err := os.WriteFile(p, []byte("FAKE charge density for testing\n"), 0o644); err != nil {
		return err
	}
	if os.Getenv("FDVIB_FAKE_PAW") != "" {
		return os.WriteFile(filepath.Join(filepath.Dir(p), "paw.txt"),
			[]byte("FAKE PAW data for testing\n"), 0o644)
	}
	return nil
}
