// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package results locates and parses the results directory metadata,
// mirroring results.cpp.
package results

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
)

// ResultMetadata is the parsed results/metadata.dat, like
// struct ResultMetadata.
type ResultMetadata struct {
	Program                 string
	ElectronicEnergyHartree float64
	Multiplicity            int
	ModeSelection           string
	SelectedAtoms           string
}

func newResultMetadata() *ResultMetadata {
	return &ResultMetadata{Program: "qe", Multiplicity: 1,
		ModeSelection: "all", SelectedAtoms: "all"}
}

func countMessage(paths []string) string {
	if len(paths) == 0 {
		return "0"
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%d (", len(paths))
	for i, p := range paths {
		if i > 0 {
			b.WriteString(", ")
		}
		b.WriteString(filepath.Base(p))
	}
	b.WriteString(")")
	return b.String()
}

// ReadMetadata parses results/metadata.dat, like result_metadata().
func ReadMetadata(results string, required bool) (*ResultMetadata, error) {
	metadata := newResultMetadata()
	path := filepath.Join(results, "metadata.dat")
	if _, err := os.Stat(path); err != nil {
		if required {
			return nil, fmt.Errorf("metadata.dat is required in %s", results)
		}
		return metadata, nil
	}
	c, err := config.Load(path)
	if err != nil {
		return nil, err
	}
	if err := c.RequireOnly(map[string]bool{
		"program": true, "electronic_energy_hartree": true,
		"multiplicity": true, "mode_selection": true, "selected_atoms": true,
	}, "metadata.dat"); err != nil {
		return nil, err
	}
	metadata.Program = config.Lower(c.Get("program", metadata.Program))
	if metadata.Program != "qe" {
		return nil, fmt.Errorf("unsupported program in metadata.dat: %s", metadata.Program)
	}
	metadata.ElectronicEnergyHartree, err = c.Real("electronic_energy_hartree", metadata.ElectronicEnergyHartree)
	if err != nil {
		return nil, err
	}
	metadata.Multiplicity, err = c.Integer("multiplicity", metadata.Multiplicity)
	if err != nil {
		return nil, err
	}
	metadata.ModeSelection = config.Lower(c.Get("mode_selection", metadata.ModeSelection))
	metadata.SelectedAtoms = config.Lower(c.Get("selected_atoms", metadata.SelectedAtoms))
	if metadata.ModeSelection != "all" && metadata.ModeSelection != "gas" &&
		metadata.ModeSelection != "local" {
		return nil, fmt.Errorf("mode_selection must be all, gas, or local in metadata.dat")
	}
	if math.IsInf(metadata.ElectronicEnergyHartree, 0) ||
		math.IsNaN(metadata.ElectronicEnergyHartree) || metadata.Multiplicity <= 0 {
		return nil, fmt.Errorf("invalid metadata.dat")
	}
	return metadata, nil
}

// ResultFiles finds the single .dynG and .freq.out in the results directory,
// like result_files().
func ResultFiles(results, context string) (dyn, freq string, err error) {
	info, err := os.Stat(results)
	if err != nil || !info.IsDir() {
		return "", "", fmt.Errorf("not a result directory: %s", results)
	}
	var dyns, freqs []string
	entries, err := os.ReadDir(results)
	if err != nil {
		return "", "", fmt.Errorf("not a result directory: %s", results)
	}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()
		if strings.HasSuffix(name, ".dynG") {
			dyns = append(dyns, filepath.Join(results, name))
		}
		if strings.HasSuffix(name, ".freq.out") {
			freqs = append(freqs, filepath.Join(results, name))
		}
	}
	if len(dyns) != 1 || len(freqs) != 1 {
		return "", "", fmt.Errorf("%s requires exactly one .dynG and one .freq.out in %s; found .dynG=%s, .freq.out=%s",
			context, results, countMessage(dyns), countMessage(freqs))
	}
	return dyns[0], freqs[0], nil
}
