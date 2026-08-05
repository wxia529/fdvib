// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package qeinput

import (
	"github.com/wxia529/fdvib/internal/fixtures"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeScf(t *testing.T, dir, content string) string {
	t.Helper()
	p := filepath.Join(dir, "scf.in")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestParseQeInput(t *testing.T) {
	q, err := ParseQeInput(writeScf(t, t.TempDir(), fixtures.WaterScf))
	if err != nil {
		t.Fatal(err)
	}
	if q.Nat != 3 || q.Ntyp != 2 {
		t.Errorf("nat/ntyp = %d/%d", q.Nat, q.Ntyp)
	}
	if q.Prefix != "h2o" {
		t.Errorf("prefix = %q", q.Prefix)
	}
	if q.PositionsUnit != "angstrom" || q.CellUnit != "angstrom" {
		t.Errorf("units = %q/%q", q.PositionsUnit, q.CellUnit)
	}
	if q.AlatAngstrom != 10.0 {
		t.Errorf("alat = %v", q.AlatAngstrom)
	}
	if len(q.Species) != 2 || q.Species[0].Symbol != "O" || q.Species[1].Mass != 1.00794 {
		t.Errorf("species: %+v", q.Species)
	}
	if len(q.Atoms) != 3 || q.Atoms[1].Type != 2 || q.Atoms[1].R[0] != 0.757 {
		t.Errorf("atoms: %+v", q.Atoms)
	}
	if q.Cell[2][2] != 10.0 {
		t.Errorf("cell: %+v", q.Cell)
	}
}

func TestParseQeInputUnits(t *testing.T) {
	// bohr positions, alat cell with celldm(1)
	content := strings.ReplaceAll(fixtures.WaterScf,
		"ATOMIC_POSITIONS angstrom\n  O  0.0  0.0  0.0\n",
		"ATOMIC_POSITIONS bohr\n  O  0.0  0.0  0.0\n")
	content = strings.ReplaceAll(content,
		"CELL_PARAMETERS angstrom\n", "CELL_PARAMETERS alat\n")
	content = strings.ReplaceAll(content, "A = 10.0", "celldm(1) = 10.0")
	q, err := ParseQeInput(writeScf(t, t.TempDir(), content))
	if err != nil {
		t.Fatal(err)
	}
	if q.AlatAngstrom != 10.0*0.529177210544 {
		t.Errorf("alat = %v", q.AlatAngstrom)
	}
	// CELL_PARAMETERS alat: vectors are multiples of alat.
	if q.Cell[0][0] != 10.0*q.AlatAngstrom {
		t.Errorf("alat cell not scaled: %+v", q.Cell)
	}
	// bohr positions: H at 0.757 bohr = 0.757 * BOHR_TO_ANG angstrom.
	if math.Abs(q.Atoms[1].R[0]-0.757*0.529177210544) > 1e-15 {
		t.Errorf("bohr position not scaled: %+v", q.Atoms[1].R)
	}
}

func TestParseQeInputRejections(t *testing.T) {
	cases := map[string]string{
		"ibrav":      strings.ReplaceAll(fixtures.WaterScf, "ibrav = 0", "ibrav = 1"),
		"not scf":    strings.ReplaceAll(fixtures.WaterScf, "calculation = 'scf'", "calculation = 'relax'"),
		"no tprnfor": strings.ReplaceAll(fixtures.WaterScf, "tprnfor = .true.", ""),
		"startingpot": strings.ReplaceAll(fixtures.WaterScf, "tprnfor = .true.",
			"tprnfor = .true.\n  startingpot = 'file'"),
		"no unit":         strings.ReplaceAll(fixtures.WaterScf, "ATOMIC_POSITIONS angstrom", "ATOMIC_POSITIONS"),
		"bad unit":        strings.ReplaceAll(fixtures.WaterScf, "ATOMIC_POSITIONS angstrom", "ATOMIC_POSITIONS car"),
		"crystal_sg":      strings.ReplaceAll(fixtures.WaterScf, "ATOMIC_POSITIONS angstrom", "ATOMIC_POSITIONS crystal_sg"),
		"unknown species": strings.ReplaceAll(fixtures.WaterScf, "H  0.757", "X  0.757"),
	}
	for name, content := range cases {
		if _, err := ParseQeInput(writeScf(t, t.TempDir(), content)); err == nil {
			t.Errorf("case %s should fail", name)
		}
	}
}

func TestReferenceInput(t *testing.T) {
	q, err := ParseQeInput(writeScf(t, t.TempDir(), fixtures.WaterScf))
	if err != nil {
		t.Fatal(err)
	}
	got, err := ReferenceInput(q, "./out", "/case", "/run/init_scf_001")
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{
		"outdir = './out',\n",
		"pseudo_dir = '../../pseudo',\n",
		"  H  0.757  0.586  0.0\n",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("missing %q in:\n%s", want, got)
		}
	}
}

func TestReferenceInputHomePseudo(t *testing.T) {
	content := strings.ReplaceAll(fixtures.WaterScf, "pseudo_dir = '../pseudo'", "pseudo_dir = '~/lib/pseudo'")
	q, err := ParseQeInput(writeScf(t, t.TempDir(), content))
	if err != nil {
		t.Fatal(err)
	}
	t.Setenv("HOME", "/home/test")
	got, err := ReferenceInput(q, "./out", "/case", "/run/x")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(got, "pseudo_dir = '/home/test/lib/pseudo',\n") {
		t.Errorf("home expansion failed:\n%s", got)
	}
}

func TestDisplacedInput(t *testing.T) {
	q, err := ParseQeInput(writeScf(t, t.TempDir(), fixtures.WaterScf))
	if err != nil {
		t.Fatal(err)
	}
	// atom is 0-based; displace atom 1 (O) along x by +0.01 angstrom.
	got, err := DisplacedInput(q, 0, 0, 0.01, "./out", "/case", "/run/disp_0001_x_p_001")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(got, "O          0.0100000000") {
		t.Errorf("displacement missing:\n%s", got)
	}
	for _, want := range []string{
		"disk_io = 'minimal',\n",
		"startingpot = 'file',\n",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("missing %q in:\n%s", want, got)
		}
	}
	// Unperturbed atoms keep their original card text.
	if !strings.Contains(got, "  H  0.757  0.586  0.0\n") {
		t.Errorf("original card text changed:\n%s", got)
	}
}

func TestDisplacedInputNoElectrons(t *testing.T) {
	content := strings.ReplaceAll(fixtures.WaterScf, "&ELECTRONS\n  conv_thr = 1.0D-8\n/\n", "")
	q, err := ParseQeInput(writeScf(t, t.TempDir(), content))
	if err != nil {
		t.Fatal(err)
	}
	got, err := DisplacedInput(q, 0, 2, -0.01, "./out", "/case", "/run/d")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(got, "&ELECTRONS\n  startingpot = 'file',\n/\n") {
		t.Errorf("synthetic &ELECTRONS missing:\n%s", got)
	}
	if !strings.Contains(got, "-0.0100000000") {
		t.Errorf("displacement wrong:\n%s", got)
	}
}

func TestDisplacedInputCrystal(t *testing.T) {
	content := strings.ReplaceAll(fixtures.WaterScf, "ATOMIC_POSITIONS angstrom", "ATOMIC_POSITIONS crystal")
	content = strings.ReplaceAll(content, "  O  0.0  0.0  0.0\n", "  O  0.5  0.5  0.5\n")
	content = strings.ReplaceAll(content, "  H  0.757  0.586  0.0\n", "  H  0.25  0.5  0.5\n")
	content = strings.ReplaceAll(content, "  H  -0.757  0.586  0.0\n", "  H  0.75  0.5  0.5\n")
	q, err := ParseQeInput(writeScf(t, t.TempDir(), content))
	if err != nil {
		t.Fatal(err)
	}
	// Displace O by 1 angstrom along x: fractional delta = 0.1.
	got, err := DisplacedInput(q, 0, 0, 1.0, "./out", "/case", "/run/d")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(got, "0.6000000000") {
		t.Errorf("crystal displacement wrong:\n%s", got)
	}
}

func TestParseQeInputIstreamSemantics(t *testing.T) {
	// C++ istream >> double uses strtod prefix parsing: "15.9994D0"
	// parses as 15.9994 and a D exponent in an atom coordinate makes the
	// whole line fail, exactly like the C++ implementation.
	content := strings.ReplaceAll(fixtures.WaterScf, "O  15.9994  O.pbe.UPF", "O  15.9994D0  O.pbe.UPF")
	q, err := ParseQeInput(writeScf(t, t.TempDir(), content))
	if err != nil {
		t.Fatal(err)
	}
	if q.Species[0].Mass != 15.9994 {
		t.Errorf("mass = %v, want 15.9994", q.Species[0].Mass)
	}
	content = strings.ReplaceAll(fixtures.WaterScf, "H  0.757  0.586  0.0", "H  0.757D0  0.586  0.0")
	if _, err := ParseQeInput(writeScf(t, t.TempDir(), content)); err == nil {
		t.Error("D exponent in atom coordinates must fail like C++")
	}
	content = strings.ReplaceAll(fixtures.WaterScf, "10.0  0.0  0.0", "10.0D0  0.0  0.0")
	if _, err := ParseQeInput(writeScf(t, t.TempDir(), content)); err == nil {
		t.Error("D exponent in cell parameters must fail like C++")
	}
}

func TestParseQeInputMissingCards(t *testing.T) {
	// Removing ATOMIC_POSITIONS must fail with the block-requirement error
	// (C++ pos_start stays -1), not a generic parse failure.
	content := strings.ReplaceAll(fixtures.WaterScf,
		"ATOMIC_POSITIONS angstrom\n  O  0.0  0.0  0.0\n  H  0.757  0.586  0.0\n  H  -0.757  0.586  0.0\n", "")
	_, err := ParseQeInput(writeScf(t, t.TempDir(), content))
	if err == nil || !strings.Contains(err.Error(), "require ATOMIC_SPECIES, ATOMIC_POSITIONS, and CELL_PARAMETERS") {
		t.Errorf("missing positions card: %v", err)
	}
	content = strings.ReplaceAll(fixtures.WaterScf, "CELL_PARAMETERS angstrom\n  10.0  0.0  0.0\n  0.0  10.0  0.0\n  0.0  0.0  10.0\n", "")
	_, err = ParseQeInput(writeScf(t, t.TempDir(), content))
	if err == nil || !strings.Contains(err.Error(), "require ATOMIC_SPECIES, ATOMIC_POSITIONS, and CELL_PARAMETERS") {
		t.Errorf("missing cell card: %v", err)
	}
}
