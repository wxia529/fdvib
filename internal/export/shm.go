// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package export

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/elements"
	modeselect "github.com/wxia529/fdvib/internal/modes"
	"github.com/wxia529/fdvib/internal/process"
	"github.com/wxia529/fdvib/internal/qevibration"
	"github.com/wxia529/fdvib/internal/results"
	"github.com/wxia529/fdvib/internal/units"
)

// shellDisplayPath renders a path for shell use, quoting when it contains
// characters outside a safe set, like shell_display_path().
func shellDisplayPath(path string) string {
	shown := config.DisplayPath(path)
	const safe = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_@%+=:,./-"
	unsafe := false
	for _, c := range shown {
		if !strings.ContainsRune(safe, c) {
			unsafe = true
			break
		}
	}
	if shown == "" || unsafe {
		return process.ShellQuote(shown)
	}
	return shown
}

// validateShm re-parses a generated .shm file, like validate_shm().
func validateShm(path string) error {
	text, err := config.ReadText(path)
	if err != nil {
		return err
	}
	lines := config.SplitLines(text)
	for _, line := range lines {
		if line == "" {
			return fmt.Errorf("generated SHM contains an empty line")
		}
	}
	if len(lines) == 0 || text[len(text)-1] != '\n' {
		return fmt.Errorf("generated SHM must end with a newline")
	}
	names := []string{"*E", "*wavenum", "*atoms", "*elevel"}
	tags := make([]int, 4)
	for ti, name := range names {
		var found []int
		for i, line := range lines {
			if strings.HasPrefix(line, name) {
				found = append(found, i)
			}
		}
		if len(found) != 1 {
			return fmt.Errorf("generated SHM must contain exactly one %s tag", name)
		}
		tags[ti] = found[0]
	}
	if !(tags[0]+1 < tags[1] && tags[1] < tags[2] && tags[2]+1 < tags[3]) {
		return fmt.Errorf("generated SHM sections are missing data or out of order")
	}
	if tags[1] != tags[0]+2 {
		return fmt.Errorf("generated SHM *E section must contain one value")
	}
	{
		fields := strings.Fields(lines[tags[0]+1])
		value, err := config.ParseNumber(fields[0])
		if len(fields) != 1 || err != nil {
			return fmt.Errorf("invalid SHM electronic energy")
		}
		_ = value
	}
	for i := tags[1] + 1; i < tags[2]; i++ {
		fields := strings.Fields(lines[i])
		if _, err := config.ParseNumber(fields[0]); len(fields) != 1 || err != nil {
			return fmt.Errorf("invalid SHM wavenumber")
		}
	}
	for i := tags[2] + 1; i < tags[3]; i++ {
		fields := strings.Fields(lines[i])
		if len(fields) != 5 {
			return fmt.Errorf("invalid SHM atom row")
		}
		mass, err := config.ParseNumber(fields[1])
		if err != nil {
			return fmt.Errorf("invalid SHM atom row")
		}
		coords := make([]float64, 3)
		for k := 0; k < 3; k++ {
			x, err := config.ParseNumber(fields[2+k])
			if err != nil {
				return fmt.Errorf("invalid SHM atom row")
			}
			coords[k] = x
		}
		if !elements.IsStandardElement(fields[0]) || !(mass > 0.0) ||
			math.IsInf(mass, 0) || math.IsNaN(mass) ||
			math.IsInf(coords[0], 0) || math.IsNaN(coords[0]) ||
			math.IsInf(coords[1], 0) || math.IsNaN(coords[1]) ||
			math.IsInf(coords[2], 0) || math.IsNaN(coords[2]) {
			return fmt.Errorf("invalid SHM atom row")
		}
	}
	hasGround := false
	for i := tags[3] + 1; i < len(lines); i++ {
		fields := strings.Fields(lines[i])
		if len(fields) != 2 {
			return fmt.Errorf("invalid SHM electronic-level row")
		}
		energy, err := config.ParseNumber(fields[0])
		if err != nil {
			return fmt.Errorf("invalid SHM electronic-level row")
		}
		degeneracy, err := strconv.Atoi(fields[1])
		if err != nil {
			return fmt.Errorf("invalid SHM electronic-level row")
		}
		if math.IsInf(energy, 0) || math.IsNaN(energy) || energy < 0.0 || degeneracy <= 0 {
			return fmt.Errorf("invalid SHM electronic-level row")
		}
		if energy == 0.0 {
			hasGround = true
		}
	}
	if tags[3]+1 == len(lines) || !hasGround {
		return fmt.Errorf("generated SHM requires a ground electronic level")
	}
	return nil
}

// Shm runs the `shm` command, like shm() in shm.cpp.
func Shm(resultsDir string) error {
	dyn, freq, err := results.ResultFiles(resultsDir, "SHM export")
	if err != nil {
		return err
	}
	geometry, err := qevibration.ReadQeDynGeometry(dyn)
	if err != nil {
		return err
	}
	modes, err := qevibration.ReadQeDynmatModes(freq, len(geometry.Masses))
	if err != nil {
		return err
	}
	metadata, err := results.ReadMetadata(resultsDir, false)
	if err != nil {
		return err
	}
	energy := metadata.ElectronicEnergyHartree
	multiplicity := metadata.Multiplicity
	selected, err := modeselect.SelectShmModes(metadata, geometry, modes)
	if err != nil {
		return err
	}

	smallestRetained := 0.0
	if len(selected.Indices) > 0 {
		smallestRetained = math.Inf(1)
		for _, index := range selected.Indices {
			if math.Abs(modes[index].Freq) < smallestRetained {
				smallestRetained = math.Abs(modes[index].Freq)
			}
		}
	}
	outputSymbols := make([]string, len(geometry.Symbols))
	for i := range geometry.Symbols {
		sym, err := elements.StandardElementSymbol(geometry.Symbols[i], "SHM export")
		if err != nil {
			return err
		}
		outputSymbols[i] = sym
		if !(geometry.Masses[i] > 0.0) {
			return fmt.Errorf("shm export requires positive atomic masses")
		}
	}

	var b strings.Builder
	b.WriteString("*E  //Electronic energy (a.u.)\n")
	b.WriteString(config.FormatSciUpper(energy, 15))
	b.WriteString("\n*wavenum  //Wavenumbers (cm-1). Negative value means imaginary frequency\n")
	for _, index := range selected.Indices {
		b.WriteString(config.FormatFixed(modes[index].Freq, 10))
		b.WriteString("\n")
	}
	b.WriteString("*atoms  //Information of all atoms: Name, mass (amu), X, Y, Z (Angstrom)\n")
	for i := range geometry.Symbols {
		b.WriteString(outputSymbols[i] + " " +
			config.FormatFixed(geometry.Masses[i], 10) + " " +
			config.FormatFixed(geometry.RBohr[i][0]*units.BohrToAng, 10) + " " +
			config.FormatFixed(geometry.RBohr[i][1]*units.BohrToAng, 10) + " " +
			config.FormatFixed(geometry.RBohr[i][2]*units.BohrToAng, 10) + "\n")
	}
	b.WriteString("*elevel  //Energy (eV) and degeneracy of electronic energy levels\n")
	fmt.Fprintf(&b, "0.000000 %d\n", multiplicity)

	destination := filepath.Join(resultsDir, strings.TrimSuffix(filepath.Base(dyn), ".dynG")+".shm")
	temporary := destination + ".tmp"
	if _, err := os.Stat(temporary); err == nil {
		return fmt.Errorf("stale SHM temporary file exists: %s", temporary)
	}
	if err := config.WriteText(temporary, b.String()); err != nil {
		return err
	}
	if err := validateShm(temporary); err != nil {
		return err
	}
	if err := os.Rename(temporary, destination); err != nil {
		return err
	}
	largestRemoved := selected.LargestRemoved
	fmt.Printf("Wrote %s\n", config.DisplayPath(destination))
	fmt.Printf("SHM mode selection: %s, retained %d, removed %d, largest removed |frequency|=%s, smallest retained |frequency|=%s cm^-1\n",
		selected.Classification, len(selected.Indices), len(modes)-len(selected.Indices),
		config.FormatGeneral(largestRemoved, 6), config.FormatGeneral(smallestRetained, 6))
	if metadata.ModeSelection == "local" {
		fmt.Printf("Run Shermo with: Shermo %s -imode 1 -PGlabel C1\n", shellDisplayPath(destination))
	} else {
		fmt.Printf("Run Shermo with: Shermo %s\n", shellDisplayPath(destination))
	}
	return nil
}
