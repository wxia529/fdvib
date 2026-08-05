// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package qeoutput

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/wxia529/fdvib/internal/config"
)

const pwOutOK = `     Program PWSCF v.7.2

     Self-consistency achieved

     Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000     0.00000000     0.00000000
     atom    2 type  2   force =     0.00000123D+00    -0.00000234     0.00000345
     atom    3 type  2   force =    -0.00000123     0.00000234    -0.00000345

     Total force =     0.00000417     Total stress =      0.00000000

!    total energy              =     -17.12345678901234 Ry
     JOB DONE.
`

func write(t *testing.T, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "out.txt")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestParseForces(t *testing.T) {
	f, err := ParseForces(write(t, pwOutOK), 3)
	if err != nil {
		t.Fatal(err)
	}
	if f[1][0] != 0.00000123 || f[2][2] != -0.00000345 || f[0][0] != 0.0 {
		t.Errorf("forces: %+v", f)
	}
}

func TestParseForcesLastBlock(t *testing.T) {
	// Two force blocks: the last one must win.
	text := pwOutOK + "\n     Forces acting on atoms (cartesian axes, Ry/au):\n" +
		"     atom    1 type  1   force =     1.00000000     2.00000000     3.00000000\n" +
		"     atom    2 type  2   force =     4.00000000     5.00000000     6.00000000\n" +
		"     atom    3 type  2   force =     7.00000000     8.00000000     9.00000000\n" +
		"     JOB DONE.\n"
	f, err := ParseForces(write(t, text), 3)
	if err != nil {
		t.Fatal(err)
	}
	if f[0][0] != 1.0 || f[2][2] != 9.0 {
		t.Errorf("forces: %+v", f)
	}
}

func TestValidateQeOutput(t *testing.T) {
	if err := ValidateQeOutput(write(t, pwOutOK)); err != nil {
		t.Errorf("valid output rejected: %v", err)
	}
	cases := map[string]string{
		"no JOB DONE":      strings.ReplaceAll(pwOutOK, "JOB DONE.", "killed"),
		"not converged":    strings.ReplaceAll(pwOutOK, "JOB DONE.", "convergence NOT achieved"),
		"error in routine": strings.ReplaceAll(pwOutOK, "JOB DONE.", "Error in routine c_bands"),
	}
	for name, content := range cases {
		if err := ValidateQeOutput(write(t, content)); err == nil {
			t.Errorf("case %s should fail", name)
		}
	}
}

func TestParseForcesErrors(t *testing.T) {
	cases := map[string]string{
		"no block": strings.ReplaceAll(pwOutOK, "Forces acting on atoms", "Forces acting"),
		"duplicate": strings.ReplaceAll(pwOutOK, "atom    1 type  1   force =     0.00000000     0.00000000     0.00000000",
			"atom    1 type  1   force =     0.00000000     0.00000000     0.00000000\n     atom    1 type  1   force =     0.00000000     0.00000000     0.00000000"),
		"incomplete": strings.ReplaceAll(pwOutOK,
			"     atom    3 type  2   force =    -0.00000123     0.00000234    -0.00000345", ""),
	}
	for name, content := range cases {
		if _, err := ParseForces(write(t, content), 3); err == nil {
			t.Errorf("case %s should fail", name)
		}
	}
}

func TestReadTotalEnergyHartree(t *testing.T) {
	e, err := ReadTotalEnergyHartree(write(t, pwOutOK))
	if err != nil {
		t.Fatal(err)
	}
	if e != -17.12345678901234/2.0 {
		t.Errorf("energy = %v", e)
	}
	// Last line wins.
	text := pwOutOK + "\n!    total energy              =     -99.0 Ry\n     JOB DONE.\n"
	e, err = ReadTotalEnergyHartree(write(t, text))
	if err != nil {
		t.Fatal(err)
	}
	if e != -49.5 {
		t.Errorf("energy = %v", e)
	}
}

func TestForcesRoundTrip(t *testing.T) {
	vecs := []config.Vec3{{1.5, -2.5e-4, 0}, {0, 0, 3.25}, {-1e-15, 2.5e15, -0.0}}
	dir := t.TempDir()
	p := filepath.Join(dir, "forces.dat")
	if err := WriteForces(p, vecs, "/scratch/pw.out"); err != nil {
		t.Fatal(err)
	}
	got, err := ReadForces(p, 3)
	if err != nil {
		t.Fatal(err)
	}
	if got[0][1] != vecs[0][1] || got[2][2] != vecs[2][2] {
		t.Errorf("round trip: %+v", got)
	}
	// Comment and index parsing.
	data, _ := os.ReadFile(p)
	text := string(data)
	for _, want := range []string{"# source: pw.out", "# nat: 3", "# units: Ry/Bohr", "# atom fx fy fz", "1 "} {
		if !strings.Contains(text, want) {
			t.Errorf("missing %q in:\n%s", want, text)
		}
	}
}

func TestReadForcesErrors(t *testing.T) {
	p := filepath.Join(t.TempDir(), "forces.dat")
	os.WriteFile(p, []byte("# only comment\n"), 0o644)
	if _, err := ReadForces(p, 1); err == nil {
		t.Error("empty forces.dat should fail")
	}
	os.WriteFile(p, []byte("1 1.0 2.0\n"), 0o644)
	if _, err := ReadForces(p, 1); err == nil {
		t.Error("short row should fail")
	}
}
