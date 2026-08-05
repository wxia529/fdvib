// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package qeoutput validates pw.x stdout and extracts forces and the total
// energy, mirroring qe_output.cpp.
package qeoutput

import (
	"fmt"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
)

var (
	convergenceRe  = regexp.MustCompile(`(?i)convergence\s+NOT\s+achieved`)
	errorRoutineRe = regexp.MustCompile(`(?i)Error\s+in\s+routine`)
	forceLineRe    = regexp.MustCompile(`(?i)^\s*atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)`)
	totalEnergyRe  = regexp.MustCompile(`(?i)!\s+total\s+energy\s*=\s*([-+0-9.EeDd]+)\s+Ry\b`)
)

// ValidateQeOutput checks a pw.x stdout file for completion and errors, like
// validate_qe_output().
func ValidateQeOutput(output string) error {
	text, err := config.ReadText(output)
	if err != nil {
		return err
	}
	if !strings.Contains(text, "JOB DONE") {
		return fmt.Errorf("missing 'JOB DONE' marker in %s", output)
	}
	if convergenceRe.MatchString(text) {
		return fmt.Errorf("scf convergence was not achieved in %s", output)
	}
	if errorRoutineRe.MatchString(text) {
		return fmt.Errorf("qe reported 'Error in routine' in %s", output)
	}
	return nil
}

// ParseForces extracts the last force block from a pw.x stdout file, like
// parse_forces().
func ParseForces(output string, nat int) ([]config.Vec3, error) {
	if err := ValidateQeOutput(output); err != nil {
		return nil, err
	}
	text, err := config.ReadText(output)
	if err != nil {
		return nil, err
	}
	pos := strings.LastIndex(text, "Forces acting on atoms")
	if pos < 0 {
		return nil, fmt.Errorf("force block not found in %s", output)
	}
	f := make([]config.Vec3, nat)
	got := make([]bool, nat)
	count := 0
	started := false
	for _, line := range strings.Split(text[pos:], "\n") {
		m := forceLineRe.FindStringSubmatch(line)
		if m != nil {
			started = true
			i, _ := strconv.Atoi(m[1])
			i--
			if i < 0 || i >= nat || got[i] {
				return nil, fmt.Errorf("invalid or duplicate atom in force block %s", output)
			}
			for k := 0; k < 3; k++ {
				x, err := config.Number(m[2+k])
				if err != nil {
					return nil, err
				}
				f[i][k] = x
			}
			got[i] = true
			count++
			if count == nat {
				return f, nil
			}
		} else if started {
			break
		}
	}
	return nil, fmt.Errorf("incomplete force block in %s", output)
}

// WriteForces writes a forces.dat cache file with 15-decimal scientific
// formatting, like write_forces().
func WriteForces(p string, f []config.Vec3, source string) error {
	var b strings.Builder
	b.WriteString("# source: " + filepath.Base(source) + "\n")
	fmt.Fprintf(&b, "# nat: %d\n", len(f))
	b.WriteString("# units: Ry/Bohr\n# atom fx fy fz\n")
	for i, v := range f {
		fmt.Fprintf(&b, "%d %s %s %s\n", i+1,
			config.FormatSci(v[0], 15), config.FormatSci(v[1], 15), config.FormatSci(v[2], 15))
	}
	return config.WriteText(p, b.String())
}

// ReadForces reads a forces.dat cache file, like read_forces().
func ReadForces(p string, nat int) ([]config.Vec3, error) {
	text, err := config.ReadText(p)
	if err != nil {
		return nil, err
	}
	f := make([]config.Vec3, nat)
	got := make([]bool, nat)
	for _, line := range strings.Split(text, "\n") {
		line = config.Trim(line)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 4 {
			return nil, fmt.Errorf("malformed forces.dat %s", p)
		}
		i, err := strconv.Atoi(fields[0])
		if err != nil || i < 1 || i > nat {
			return nil, fmt.Errorf("malformed forces.dat %s", p)
		}
		for k := 0; k < 3; k++ {
			// C++ reads with istream >> double (strtod prefix semantics);
			// an unconsumed remainder fails the following read.
			x, rest, ok := config.IstreamDouble(fields[1+k])
			if !ok || rest != "" {
				return nil, fmt.Errorf("malformed forces.dat %s", p)
			}
			f[i-1][k] = x
		}
		got[i-1] = true
	}
	for i := 0; i < nat; i++ {
		if !got[i] {
			return nil, fmt.Errorf("incomplete forces.dat %s", p)
		}
	}
	return f, nil
}

// ReadTotalEnergyHartree extracts the last converged QE total energy in
// Hartree, like read_total_energy_hartree().
func ReadTotalEnergyHartree(output string) (float64, error) {
	text, err := config.ReadText(output)
	if err != nil {
		return 0, err
	}
	found := false
	var energyRy float64
	for _, m := range totalEnergyRe.FindAllStringSubmatch(text, -1) {
		x, err := config.Number(m[1])
		if err != nil {
			return 0, err
		}
		energyRy = x
		found = true
	}
	if !found {
		return 0, fmt.Errorf("cannot find converged qe total energy in %s", output)
	}
	return energyRy / 2.0, nil
}
