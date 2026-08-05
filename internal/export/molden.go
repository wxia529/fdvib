// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package export writes the Molden (.mol), Shermo (.shm), and thermochemistry
// (thermo.dat) output files with the exact C++ iostream layouts, mirroring
// molden.cpp, shm.cpp, and thermo.cpp.
package export

import (
	"fmt"
	"path/filepath"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/elements"
	"github.com/wxia529/fdvib/internal/modes"
	"github.com/wxia529/fdvib/internal/qevibration"
	"github.com/wxia529/fdvib/internal/results"
	"github.com/wxia529/fdvib/internal/units"
)

// writeMolden renders a Molden-format file, like write_molden().
func writeMolden(path string, geometry *qevibration.DynGeometry,
	modes []qevibration.Mode, selected *modes.ModeSelection) error {
	symbols := make([]string, len(geometry.Symbols))
	for i, label := range geometry.Symbols {
		sym, err := elements.StandardElementSymbol(label, "Molden export")
		if err != nil {
			return err
		}
		symbols[i] = sym
	}

	var b strings.Builder
	b.WriteString("[Molden Format]\n")
	b.WriteString("[Cell]\n")
	for _, v := range geometry.CellBohr {
		for k := 0; k < 3; k++ {
			b.WriteString(config.FixedField(v[k]*units.BohrToAng, 8, 15))
		}
		b.WriteString("\n")
	}
	b.WriteString("[Atoms] AU\n")
	for i := range symbols {
		z, err := elements.AtomicNumberFromLabel(symbols[i], "Molden export")
		if err != nil {
			return err
		}
		b.WriteString(fmt.Sprintf("%6s", symbols[i]))
		b.WriteString(fmt.Sprintf("%8d", i+1))
		b.WriteString(fmt.Sprintf("%8d", z))
		for k := 0; k < 3; k++ {
			b.WriteString(config.FixedField(geometry.RBohr[i][k], 8, 18))
		}
		b.WriteString("\n")
	}
	b.WriteString("[FREQ]\n")
	for _, index := range selected.Indices {
		b.WriteString(config.FormatFixed(modes[index].Freq, 8))
		b.WriteString("\n")
	}
	b.WriteString("[FR-COORD]\n")
	for i := range symbols {
		b.WriteString(fmt.Sprintf("%6s", symbols[i]))
		for k := 0; k < 3; k++ {
			b.WriteString(config.FixedField(geometry.RBohr[i][k], 8, 18))
		}
		b.WriteString("\n")
	}
	b.WriteString("[FR-NORM-COORD]\n")
	for i, index := range selected.Indices {
		fmt.Fprintf(&b, " vibration %d\n", i+1)
		for _, displacement := range modes[index].Displacement {
			for k := 0; k < 3; k++ {
				b.WriteString(config.FixedField(displacement[k], 10, 15))
			}
			b.WriteString("\n")
		}
	}
	if err := config.WriteText(path, b.String()); err != nil {
		return err
	}
	fmt.Printf("Wrote %d Molden modes to %s\n", len(selected.Indices), config.DisplayPath(path))
	return nil
}

// Modes runs the `modes` command, like modes() in molden.cpp.
func Modes(resultsDir string) error {
	dyn, freq, err := results.ResultFiles(resultsDir, "Normal-mode analysis")
	if err != nil {
		return err
	}
	geometry, err := qevibration.ReadQeDynGeometry(dyn)
	if err != nil {
		return err
	}
	parsed, err := qevibration.ReadQeDynmatModes(freq, len(geometry.Masses))
	if err != nil {
		return err
	}
	metadata, err := results.ReadMetadata(resultsDir, false)
	if err != nil {
		return err
	}
	selected, err := modes.SelectMoldenModes(metadata, geometry, parsed)
	if err != nil {
		return err
	}
	stem := strings.TrimSuffix(filepath.Base(dyn), ".dynG")
	return writeMolden(filepath.Join(resultsDir, stem+".mol"), geometry, parsed, selected)
}
