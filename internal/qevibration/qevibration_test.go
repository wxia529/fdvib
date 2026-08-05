// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package qevibration

import (
	"github.com/wxia529/fdvib/internal/fixtures"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestReadQeDynmatModes(t *testing.T) {
	modes, err := ReadQeDynmatModes(writeFile(t, t.TempDir(), "freq.out", fixtures.FreqOut), 3)
	if err != nil {
		t.Fatal(err)
	}
	if len(modes) != 9 {
		t.Fatalf("modes = %d", len(modes))
	}
	if modes[0].Freq != 296.679370 || modes[8].Freq != 1568.590178 {
		t.Errorf("freqs: %v %v", modes[0].Freq, modes[8].Freq)
	}
	if len(modes[0].Displacement) != 3 || modes[0].Displacement[0][0] != 0.0 ||
		modes[0].Displacement[0][1] != 0.5987416 ||
		modes[0].Displacement[1][2] != -0.4690194 {
		t.Errorf("displacements: %+v", modes[0].Displacement)
	}
}

func TestReadQeDynmatModesErrors(t *testing.T) {
	dir := t.TempDir()
	// Only 8 modes.
	bad := strings.Repeat(strings.Join(strings.Split(fixtures.FreqOut, "\n")[:4], "\n")+"\n", 8)
	if _, err := ReadQeDynmatModes(writeFile(t, dir, "short.out", bad), 3); err == nil {
		t.Error("short mode count should fail")
	}
	if _, err := ReadQeDynmatModes(writeFile(t, dir, "zero.out", fixtures.FreqOut), 0); err == nil {
		t.Error("nat=0 should fail")
	}
}

func TestReadQeDynGeometry(t *testing.T) {
	g, err := ReadQeDynGeometry(writeFile(t, t.TempDir(), "system.dynG", fixtures.DynG))
	if err != nil {
		t.Fatal(err)
	}
	if len(g.Masses) != 3 || len(g.Symbols) != 3 {
		t.Fatalf("geometry size: %+v", g)
	}
	if g.Symbols[0] != "O" || g.Masses[1] != 918.428582136069/911.4442431390707 {
		t.Errorf("species: %v %v", g.Symbols, g.Masses)
	}
	if g.CellBohr[0][0] != 10.0 {
		t.Errorf("cell: %+v", g.CellBohr)
	}
	// dynG coordinates are in alat units (alat = 10 bohr here).
	if abs(g.RBohr[1][0]-1.430510998*10.0) > 1e-9 {
		t.Errorf("r: %+v", g.RBohr)
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func TestReadQeDynGeometryErrors(t *testing.T) {
	dir := t.TempDir()
	cases := map[string]string{
		"short header": "x\ny\n",
		"ibrav":        strings.Replace(fixtures.DynG, "  2  3  0", "  2  3  1", 1),
		"bad mass":     strings.Replace(fixtures.DynG, "14578.211275740000", "nan", 1),
	}
	for name, content := range cases {
		if _, err := ReadQeDynGeometry(writeFile(t, dir, name+".dynG", content)); err == nil {
			t.Errorf("case %s should fail", name)
		}
	}
}
