// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package init

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestInitializeInputLocal(t *testing.T) {
	dir := t.TempDir()
	if err := InitializeInput("local", dir); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(filepath.Join(dir, "fdvib.in"))
	if err != nil {
		t.Fatal(err)
	}
	text := string(data)
	for _, want := range []string{
		"scf_input = scf.in\n", "outdir = fdvib\n", "system_type = local\n",
		"selected_atoms = 1,2,3\n", "displacement_angstrom = 0.01\n",
		"pw_command = pw.x\n", "prefix = system\n", "run_dynmat = true\n",
		"dynmat_command = dynmat.x\n",
	} {
		if !strings.Contains(text, want) {
			t.Errorf("missing %q in template:\n%s", want, text)
		}
	}
}

func TestInitializeInputGas(t *testing.T) {
	dir := t.TempDir()
	if err := InitializeInput("GAS", dir); err != nil {
		t.Fatal(err)
	}
	data, _ := os.ReadFile(filepath.Join(dir, "fdvib.in"))
	text := string(data)
	for _, want := range []string{"system_type = gas\n", "selected_atoms = all\n",
		"multiplicity = 1\n", "prefix = molecule\n"} {
		if !strings.Contains(text, want) {
			t.Errorf("missing %q in template:\n%s", want, text)
		}
	}
}

func TestInitializeInputRefusesOverwrite(t *testing.T) {
	dir := t.TempDir()
	if err := InitializeInput("local", dir); err != nil {
		t.Fatal(err)
	}
	if err := InitializeInput("local", dir); err == nil {
		t.Error("overwrite should be refused")
	}
}

func TestInitializeInputBadType(t *testing.T) {
	if err := InitializeInput("solid", t.TempDir()); err == nil {
		t.Error("bad type should fail")
	}
}
