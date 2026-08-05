// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package qevibration

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

const freqOut = `   freq (    1) =       8.894417 [THz] =     296.679370 [cm-1]
(   0.0000000   0.0000000   0.5987416   0.0000000  -0.1548613   0.0000000 )
(  -0.1834352   0.0000000  -0.2030341   0.0000000  -0.4690194   0.0000000 )
(  -0.4183367   0.0000000  -0.4743769   0.0000000  -0.2601189   0.0000000 )
   freq (    2) =      10.723309 [THz] =     357.656466 [cm-1]
(   0.0000000   0.0000000   0.6937382   0.0000000   0.0000000   0.0000000 )
(  -0.0000000   0.0000000   0.3586637   0.0000000  -0.4892037   0.0000000 )
(  -0.0000000   0.0000000  -0.4739768   0.0000000  -0.5094611   0.0000000 )
   freq (    3) =      13.631979 [THz] =     454.674703 [cm-1]
(   0.0000000   0.0000000  -0.5944620   0.0000000   0.1689872   0.0000000 )
(   0.0000000   0.0000000  -0.4819955   0.0000000  -0.3501150   0.0000000 )
(   0.0000000   0.0000000  -0.4562901   0.0000000   0.3093417   0.0000000 )
   freq (    4) =      15.165869 [THz] =     505.908220 [cm-1]
(  -0.2652315   0.0000000   0.0000000   0.0000000   0.3353864   0.0000000 )
(  -0.3525340   0.0000000   0.0000000   0.0000000  -0.5689768   0.0000000 )
(  -0.5939046   0.0000000   0.0000000   0.0000000  -0.3803498   0.0000000 )
   freq (    5) =      20.278004 [THz] =     676.430641 [cm-1]
(   0.0000000   0.0000000   0.5612587   0.0000000   0.3795926   0.0000000 )
(   0.0000000   0.0000000  -0.5455749   0.0000000  -0.0025027   0.0000000 )
(   0.0000000   0.0000000  -0.3736406   0.0000000   0.2014618   0.0000000 )
   freq (    6) =      20.572647 [THz] =     686.244882 [cm-1]
(  -0.6164614   0.0000000   0.0000000   0.0000000   0.4410561   0.0000000 )
(  -0.3239216   0.0000000   0.0000000   0.0000000  -0.1516139   0.0000000 )
(  -0.1917926   0.0000000   0.0000000   0.0000000  -0.4937465   0.0000000 )
   freq (    7) =      39.339697 [THz] =    1312.365848 [cm-1]
(   0.0000000   0.0000000   0.1105471   0.0000000  -0.0000000   0.0000000 )
(   0.0000000   0.0000000  -0.4192989   0.0000000   0.6990839   0.0000000 )
(   0.0000000   0.0000000  -0.5174127   0.0000000  -0.4386184   0.0000000 )
   freq (    8) =      46.379012 [THz] =    1547.100791 [cm-1]
(   0.5472858   0.0000000   0.0000000   0.0000000  -0.0971021   0.0000000 )
(  -0.5854904   0.0000000   0.0000000   0.0000000  -0.2508146   0.0000000 )
(  -0.2144412   0.0000000   0.0000000   0.0000000  -0.3560562   0.0000000 )
   freq (    9) =      47.024167 [THz] =    1568.590178 [cm-1]
(   0.5615321   0.0000000   0.0000000   0.0000000   0.3008486   0.0000000 )
(  -0.4117948   0.0000000   0.0000000   0.0000000  -0.1530787   0.0000000 )
(   0.0785855   0.0000000   0.0000000   0.0000000  -0.5488047   0.0000000 )
`

const dynG = `Dynamical matrix file
fdvib finite-difference local/Gamma matrix
  2  3  0    10.0000000           0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
Basis vectors
    1.000000000     0.000000000     0.000000000
    0.000000000     1.000000000     0.000000000
    0.000000000     0.000000000     1.000000000
  1  'O      '        14578.211275740000
  2  'H      '          918.428582136069
    1    1     0.000000000     0.000000000     0.000000000
    2    2     1.430510998     1.107412588     0.000000000
    3    2    -1.430510998     1.107412588     0.000000000

     Dynamical  Matrix in cartesian axes

     q = (    0.000000000   0.000000000   0.000000000 ) 
`

func writeFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestReadQeDynmatModes(t *testing.T) {
	modes, err := ReadQeDynmatModes(writeFile(t, t.TempDir(), "freq.out", freqOut), 3)
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
	bad := strings.Repeat(strings.Join(strings.Split(freqOut, "\n")[:4], "\n")+"\n", 8)
	if _, err := ReadQeDynmatModes(writeFile(t, dir, "short.out", bad), 3); err == nil {
		t.Error("short mode count should fail")
	}
	if _, err := ReadQeDynmatModes(writeFile(t, dir, "zero.out", freqOut), 0); err == nil {
		t.Error("nat=0 should fail")
	}
}

func TestReadQeDynGeometry(t *testing.T) {
	g, err := ReadQeDynGeometry(writeFile(t, t.TempDir(), "system.dynG", dynG))
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
		"ibrav":        strings.Replace(dynG, "  2  3  0", "  2  3  1", 1),
		"bad mass":     strings.Replace(dynG, "14578.211275740000", "nan", 1),
	}
	for name, content := range cases {
		if _, err := ReadQeDynGeometry(writeFile(t, dir, name+".dynG", content)); err == nil {
			t.Errorf("case %s should fail", name)
		}
	}
}
