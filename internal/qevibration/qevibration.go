// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package qevibration parses dynmat.x frequency output and QE dynamical
// matrix (.dynG) files, mirroring qe_vibration.cpp.
package qevibration

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/units"
)

// Mode is one normal mode with its frequency (cm^-1, signed) and per-atom
// displacements, like struct Mode.
type Mode struct {
	Freq         float64
	Displacement []config.Vec3
}

// DynGeometry is the geometry read from a .dynG file in Bohr/amu, like
// struct DynGeometry.
type DynGeometry struct {
	CellBohr [3]config.Vec3
	Masses   []float64
	RBohr    []config.Vec3
	Symbols  []string
}

var (
	frequencyRe   = regexp.MustCompile(`(?i)freq\s*\(\s*\d+\s*\)\s*=\s*[-+0-9.EeDd]+\s*\[THz\]\s*=\s*([-+0-9.EeDd]+)\s*\[cm-1\]`)
	eigenvectorRe = regexp.MustCompile(`\(\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*\)`)
	quotedRe      = regexp.MustCompile(`'([^']*)'`)
)

// ReadQeDynmatModes parses the modes in a dynmat.x frequency output file,
// like read_qe_dynmat_modes().
func ReadQeDynmatModes(path string, nat int) ([]Mode, error) {
	if nat <= 0 {
		return nil, fmt.Errorf("invalid atom count while parsing modes in %s", path)
	}
	text, err := config.ReadText(path)
	if err != nil {
		return nil, err
	}
	var modes []Mode
	lines := strings.Split(text, "\n")
	i := 0
	for ; i < len(lines); i++ {
		m := frequencyRe.FindStringSubmatch(lines[i])
		if m == nil {
			continue
		}
		freq, err := config.Number(m[1])
		if err != nil {
			return nil, err
		}
		var mode Mode
		mode.Freq = freq
		mode.Displacement = make([]config.Vec3, 0, nat)
		for a := 0; a < nat; a++ {
			i++
			if i >= len(lines) {
				return nil, fmt.Errorf("malformed eigenvector block in %s", path)
			}
			ev := eigenvectorRe.FindStringSubmatch(lines[i])
			if ev == nil {
				return nil, fmt.Errorf("malformed eigenvector block in %s", path)
			}
			var v config.Vec3
			for k := 0; k < 3; k++ {
				x, err := config.Number(ev[1+2*k])
				if err != nil {
					return nil, err
				}
				v[k] = x
			}
			mode.Displacement = append(mode.Displacement, v)
		}
		modes = append(modes, mode)
	}
	if len(modes) != 3*nat {
		return nil, fmt.Errorf("expected %d modes, found %d in %s", 3*nat, len(modes), path)
	}
	return modes, nil
}

// ReadQeDynGeometry parses the geometry in a .dynG file, like
// read_qe_dyn_geometry().
func ReadQeDynGeometry(path string) (*DynGeometry, error) {
	text, err := config.ReadText(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(text, "\n")
	// The first three lines are the file title, description, and header.
	if len(lines) < 3 {
		return nil, fmt.Errorf("incomplete dynG header in %s", path)
	}
	header := strings.Fields(lines[2])
	if len(header) < 4 {
		return nil, fmt.Errorf("unsupported dynG header")
	}
	ntyp, err1 := strconv.Atoi(header[0])
	nat, err2 := strconv.Atoi(header[1])
	ibrav, err3 := strconv.Atoi(header[2])
	alat, err4 := config.Number(header[3])
	if err1 != nil || err2 != nil || err3 != nil || err4 != nil ||
		ibrav != 0 || ntyp <= 0 || nat <= 0 || !(alat > 0.0) {
		return nil, fmt.Errorf("unsupported dynG header")
	}
	if len(lines) < 4 {
		return nil, fmt.Errorf("incomplete dynG basis section in %s", path)
	}

	g := &DynGeometry{}
	for i := 0; i < 3; i++ {
		if 4+i >= len(lines) {
			return nil, fmt.Errorf("incomplete dynG basis vectors in %s", path)
		}
		fields := strings.Fields(lines[4+i])
		if len(fields) < 3 {
			return nil, fmt.Errorf("malformed basis vector in dynG: %s", path)
		}
		for k := 0; k < 3; k++ {
			x, err := config.Number(fields[k])
			if err != nil {
				return nil, fmt.Errorf("malformed basis vector in dynG: %s", path)
			}
			g.CellBohr[i][k] = x * alat
		}
	}
	a, b, c := g.CellBohr[0], g.CellBohr[1], g.CellBohr[2]
	volume := a[0]*(b[1]*c[2]-b[2]*c[1]) -
		a[1]*(b[0]*c[2]-b[2]*c[0]) +
		a[2]*(b[0]*c[1]-b[1]*c[0])
	if math.IsNaN(volume) || math.IsInf(volume, 0) || math.Abs(volume) < 1.0e-12 {
		return nil, fmt.Errorf("singular basis vectors in dynG: %s", path)
	}

	typeMasses := make([]float64, ntyp)
	typeSymbols := make([]string, ntyp)
	for i := 0; i < ntyp; i++ {
		lineIndex := 7 + i
		if lineIndex >= len(lines) {
			return nil, fmt.Errorf("incomplete dynG species block in %s", path)
		}
		m := quotedRe.FindStringSubmatch(lines[lineIndex])
		if m == nil {
			return nil, fmt.Errorf("malformed species symbol in dynG")
		}
		typeSymbols[i] = config.Trim(m[1])
		fields := strings.Fields(lines[lineIndex])
		if len(fields) == 0 {
			return nil, fmt.Errorf("malformed species mass in dynG")
		}
		last, err := config.Number(fields[len(fields)-1])
		if err != nil {
			return nil, err
		}
		typeMasses[i] = last / units.AmuRy
	}

	g.Masses = make([]float64, nat)
	g.RBohr = make([]config.Vec3, nat)
	g.Symbols = make([]string, nat)
	for i := 0; i < nat; i++ {
		lineIndex := 7 + ntyp + i
		if lineIndex >= len(lines) {
			return nil, fmt.Errorf("incomplete dynG atom block in %s", path)
		}
		fields := strings.Fields(lines[lineIndex])
		if len(fields) < 5 {
			return nil, fmt.Errorf("malformed atom row in dynG: %s", path)
		}
		atomIndex, err1 := strconv.Atoi(fields[0])
		typeIndex, err2 := strconv.Atoi(fields[1])
		if err1 != nil || err2 != nil || atomIndex != i+1 || typeIndex < 1 || typeIndex > ntyp {
			return nil, fmt.Errorf("malformed atom row in dynG: %s", path)
		}
		for k := 0; k < 3; k++ {
			x, err := config.Number(fields[2+k])
			if err != nil {
				return nil, fmt.Errorf("malformed atom row in dynG: %s", path)
			}
			g.RBohr[i][k] = x * alat
		}
		g.Masses[i] = typeMasses[typeIndex-1]
		g.Symbols[i] = typeSymbols[typeIndex-1]
	}
	return g, nil
}
