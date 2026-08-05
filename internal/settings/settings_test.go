// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package settings

import (
	"os"
	"path/filepath"
	"testing"
)

func writeFile(t *testing.T, path, content string) {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
}

func TestSettingsDefaults(t *testing.T) {
	dir := t.TempDir()
	cfg := filepath.Join(dir, "fdvib.in")
	writeFile(t, cfg, "scf_input = scf.in\noutdir = fdvib\nrun_dynmat = true\n")
	s, err := From(cfg, "")
	if err != nil {
		t.Fatal(err)
	}
	if s.Root != dir || s.ScfInput != filepath.Join(dir, "scf.in") ||
		s.Workdir != filepath.Join(dir, "fdvib") {
		t.Errorf("paths: %+v", s)
	}
	if s.SystemType != "local" || s.Displacement != 0.01 || s.Multiplicity != 1 ||
		s.PWCommand != "pw.x" || s.DynmatCommand != "dynmat.x" ||
		s.OutputPrefix != "system" || !s.RunDynmat || s.SelectedAll {
		t.Errorf("defaults: %+v", s)
	}
}

func TestSettingsValidation(t *testing.T) {
	dir := t.TempDir()
	cases := []string{
		"displacement_angstrom = 0\n",
		"displacement_angstrom = -1\n",
		"multiplicity = 0\n",
		"multiplicity = -2\n",
		"prefix = a/b\n",
		"prefix = .\n",
		"prefix = ..\n",
		"prefix = \n",
		"pw_command =  \n",
		"dynmat_command = \nrun_dynmat = true\n",
		"system_type = solid\n",
		"bogus_key = 1\n",
		"scf_input = a\nscf_input = b\n",
		"run_dynmat = yes\n",
	}
	for i, content := range cases {
		cfg := filepath.Join(dir, "case", string(rune('a'+i)), "fdvib.in")
		writeFile(t, cfg, content)
		if _, err := From(cfg, ""); err == nil {
			t.Errorf("case %d should fail: %q", i, content)
		}
	}
}

func TestSettingsSelectedAtoms(t *testing.T) {
	dir := t.TempDir()
	cfg := filepath.Join(dir, "fdvib.in")
	writeFile(t, cfg, "selected_atoms = 1,2;3\nrun_dynmat = true\n")
	s, err := From(cfg, "")
	if err != nil {
		t.Fatal(err)
	}
	if len(s.Selected) != 3 || s.Selected[2] != 3 {
		t.Errorf("selected = %v", s.Selected)
	}
	writeFile(t, cfg, "selected_atoms = all\nrun_dynmat = true\n")
	s, err = From(cfg, "")
	if err != nil {
		t.Fatal(err)
	}
	if !s.SelectedAll || len(s.Selected) != 0 {
		t.Errorf("all: %+v", s)
	}
	writeFile(t, cfg, "run_dynmat = true\n")
	s, err = From(cfg, "")
	if err != nil || len(s.Selected) != 0 || s.SelectedAll {
		t.Errorf("empty list: %+v %v", s, err)
	}
}

func TestJobName(t *testing.T) {
	cases := []struct {
		atom, axis, sign int
		want             string
	}{
		{1, 0, 1, "disp_0001_x_p"},
		{12, 1, -1, "disp_0012_y_m"},
		{1234, 2, 1, "disp_1234_z_p"},
	}
	for _, c := range cases {
		if got := JobName(c.atom, c.axis, c.sign); got != c.want {
			t.Errorf("JobName(%d,%d,%d) = %q, want %q", c.atom, c.axis, c.sign, got, c.want)
		}
	}
}

func TestSettingsAbsolutePaths(t *testing.T) {
	// Absolute scf_input/outdir override the fdvib.in directory, matching
	// C++ fs::path operator/ semantics.
	dir := t.TempDir()
	absScf := filepath.Join(dir, "scf.abs")
	absOut := filepath.Join(dir, "out.abs")
	if err := os.WriteFile(absScf, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg := filepath.Join(dir, "fdvib.in")
	content := "scf_input = " + absScf + "\noutdir = " + absOut + "\nrun_dynmat = true\n"
	writeFile(t, cfg, content)
	s, err := From(cfg, "")
	if err != nil {
		t.Fatal(err)
	}
	if s.ScfInput != absScf {
		t.Errorf("ScfInput = %q, want %q", s.ScfInput, absScf)
	}
	if s.Workdir != absOut {
		t.Errorf("Workdir = %q, want %q", s.Workdir, absOut)
	}
}
