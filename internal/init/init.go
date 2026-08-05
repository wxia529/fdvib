// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package init generates starter fdvib.in templates, mirroring init.cpp.
package init

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/wxia529/fdvib/internal/config"
)

func executionSettings(prefix string) string {
	return "displacement_angstrom = 0.01\n" +
		"pw_command = pw.x\n" +
		"prefix = " + prefix + "\n" +
		"run_dynmat = true\n" +
		"dynmat_command = dynmat.x\n"
}

func inputTemplate(typeName string) (string, error) {
	common := "# Change scf_input if your QE input uses another filename.\n" +
		"scf_input = scf.in\n" +
		"outdir = fdvib\n"

	if typeName == "local" {
		return "# Replace selected_atoms with one-based QE atom indices.\n" + common +
			"system_type = local\n" +
			"selected_atoms = 1,2,3\n" + executionSettings("system"), nil
	}
	if typeName == "gas" {
		return "# Set multiplicity to the molecular spin multiplicity.\n" + common +
			"system_type = gas\n" +
			"selected_atoms = all\n" +
			"multiplicity = 1\n" + executionSettings("molecule"), nil
	}
	return "", fmt.Errorf("init type must be local or gas")
}

// InitializeInput writes a starter fdvib.in, refusing to overwrite an
// existing file, like initialize_input().
func InitializeInput(typeName, directory string) error {
	contents, err := inputTemplate(config.Lower(typeName))
	if err != nil {
		return err
	}
	destination := filepath.Join(directory, "fdvib.in")
	if _, err := os.Lstat(destination); err == nil {
		return fmt.Errorf("refusing to overwrite existing %s", config.DisplayPath(destination))
	}
	if err := config.WriteText(destination, contents); err != nil {
		return err
	}
	fmt.Printf("Wrote %s\n", config.DisplayPath(destination))
	return nil
}
