// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package results

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeFiles(t *testing.T, dir string, files map[string]string) {
	t.Helper()
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	for name, content := range files {
		if err := os.WriteFile(filepath.Join(dir, name), []byte(content), 0o644); err != nil {
			t.Fatal(err)
		}
	}
}

func TestReadMetadata(t *testing.T) {
	dir := t.TempDir()
	writeFiles(t, dir, map[string]string{
		"metadata.dat": "program = qe\nelectronic_energy_hartree = -1.5e+01\nmultiplicity = 2\nmode_selection = gas\nselected_atoms = all\n",
	})
	md, err := ReadMetadata(dir, true)
	if err != nil {
		t.Fatal(err)
	}
	if md.Program != "qe" || md.ElectronicEnergyHartree != -15.0 || md.Multiplicity != 2 ||
		md.ModeSelection != "gas" || md.SelectedAtoms != "all" {
		t.Errorf("metadata: %+v", md)
	}
}

func TestReadMetadataDefaults(t *testing.T) {
	dir := t.TempDir()
	// Missing metadata.dat: defaults when not required.
	md, err := ReadMetadata(dir, false)
	if err != nil {
		t.Fatal(err)
	}
	if md.Program != "qe" || md.ElectronicEnergyHartree != 0.0 || md.Multiplicity != 1 ||
		md.ModeSelection != "all" || md.SelectedAtoms != "all" {
		t.Errorf("defaults: %+v", md)
	}
	if _, err := ReadMetadata(dir, true); err == nil {
		t.Error("required missing metadata should fail")
	}
}

func TestReadMetadataRejections(t *testing.T) {
	dir := t.TempDir()
	cases := map[string]string{
		"bad program": "program = gaussian\n",
		"bad mode":    "mode_selection = solid\n",
		"zero mult":   "multiplicity = 0\n",
		"nan energy":  "electronic_energy_hartree = nan\n",
		"unknown key": "bogus = 1\n",
	}
	for name, content := range cases {
		sub := filepath.Join(dir, name)
		writeFiles(t, sub, map[string]string{"metadata.dat": content})
		if _, err := ReadMetadata(sub, true); err == nil {
			t.Errorf("case %s should fail", name)
		}
	}
}

func TestResultFiles(t *testing.T) {
	dir := t.TempDir()
	writeFiles(t, dir, map[string]string{
		"system.dynG": "x", "system.freq.out": "y",
	})
	dyn, freq, err := ResultFiles(dir, "test")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.HasSuffix(dyn, "system.dynG") || !strings.HasSuffix(freq, "system.freq.out") {
		t.Errorf("files: %s %s", dyn, freq)
	}

	// Multiple dynG files must fail with a listing.
	writeFiles(t, dir, map[string]string{"other.dynG": "z"})
	if _, _, err := ResultFiles(dir, "test"); err == nil {
		t.Error("two dynG should fail")
	} else if !strings.Contains(err.Error(), ".dynG=2 (") {
		t.Errorf("message: %v", err)
	}

	// Missing directory.
	if _, _, err := ResultFiles(filepath.Join(dir, "nope"), "test"); err == nil {
		t.Error("missing dir should fail")
	}
}
