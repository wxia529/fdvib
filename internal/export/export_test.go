// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package export

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Shared fixtures (same as the qevibration tests).
const fixtureDynG = `Dynamical matrix file
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

const fixtureFreqOut = `   freq (    1) =       8.894417 [THz] =     296.679370 [cm-1]
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

const fixtureMetadata = `program = qe
electronic_energy_hartree = -9.8123456789012345e+00
multiplicity = 1
mode_selection = local
selected_atoms = 1,2,3
`

const fixtureDynmatIn = `&INPUT
  fildyn='system.dynG',
  filout='system.freq.out',
  asr='no',
  remove_interaction_blocks=.true.,
/
`

// buildResults creates a results directory with local-mode fixtures.
func buildResults(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "fdvib", "results")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	files := map[string]string{
		"system.dynG":     fixtureDynG,
		"system.freq.out": fixtureFreqOut,
		"metadata.dat":    fixtureMetadata,
		"dynmat.in":       fixtureDynmatIn,
	}
	for name, content := range files {
		if err := os.WriteFile(filepath.Join(dir, name), []byte(content), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	return dir
}

func TestModesCommand(t *testing.T) {
	dir := buildResults(t)
	if err := Modes(dir); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(filepath.Join(dir, "system.mol"))
	if err != nil {
		t.Fatal(err)
	}
	text := string(data)
	for _, want := range []string{
		"[Molden Format]\n", "[Cell]\n", "[Atoms] AU\n", "[FREQ]\n",
		"[FR-COORD]\n", "[FR-NORM-COORD]\n", " vibration 1\n", " vibration 9\n",
	} {
		if !strings.Contains(text, want) {
			t.Errorf("missing %q in molden output", want)
		}
	}
	// All frequencies are nonzero in this fixture, so 9 modes are kept.
	if strings.Count(text, " vibration ") != 9 {
		t.Errorf("expected 9 vibrations, got:\n%s", text)
	}
	// [Atoms] AU positions in Bohr with fixed 8 decimals.
	fields := strings.Fields(text)
	joined := strings.Join(fields, " ")
	if !strings.Contains(joined, "O 1 8 0.00000000 0.00000000 0.00000000 H 2 1 14.30510998 11.07412588 0.00000000") {
		t.Errorf("atoms rows wrong:\n%s", text)
	}
	// [Cell] in Angstrom.
	if !strings.Contains(text, "     5.29177211") {
		t.Errorf("cell in angstrom wrong:\n%s", text)
	}
}

func TestShmCommand(t *testing.T) {
	dir := buildResults(t)
	if err := Shm(dir); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(filepath.Join(dir, "system.shm"))
	if err != nil {
		t.Fatal(err)
	}
	text := string(data)
	for _, want := range []string{
		"*E  //Electronic energy (a.u.)\n",
		"-9.812345678901234E+00\n",
		"*wavenum  //Wavenumbers (cm-1). Negative value means imaginary frequency\n",
		"*atoms  //Information of all atoms: Name, mass (amu), X, Y, Z (Angstrom)\n",
		"*elevel  //Energy (eV) and degeneracy of electronic energy levels\n",
		"0.000000 1\n",
		"O 15.9946276313 0.0000000000 0.0000000000 0.0000000000\n",
	} {
		if !strings.Contains(text, want) {
			t.Errorf("missing %q in shm output:\n%s", want, text)
		}
	}
	// local selection keeps 3*3 = 9 modes (all nonzero here).
	if strings.Count(text, "\n") != 18 {
		t.Errorf("unexpected line count: %d\n%s", strings.Count(text, "\n"), text)
	}
	// No stale tmp file left behind.
	if _, err := os.Stat(filepath.Join(dir, "system.shm.tmp")); err == nil {
		t.Error("stale tmp file remains")
	}
}

func TestShmRefusesStaleTmp(t *testing.T) {
	dir := buildResults(t)
	os.WriteFile(filepath.Join(dir, "system.shm.tmp"), []byte("x"), 0o644)
	if err := Shm(dir); err == nil {
		t.Error("stale tmp should fail")
	}
}

func TestThermoCommandLocal(t *testing.T) {
	dir := buildResults(t)
	thermoIn := filepath.Join(filepath.Dir(filepath.Dir(dir)), "thermo.in")
	content := "model = local_harmonic\ntemperature_k = 298.15\n" +
		"low_frequency_model = frequency_floor\nfrequency_floor_cm1 = 100.0\n" +
		"zero_tolerance_cm1 = 1.0\n"
	if err := os.WriteFile(thermoIn, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := Thermo(dir, thermoIn); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(filepath.Join(dir, "thermo.dat"))
	if err != nil {
		t.Fatal(err)
	}
	text := string(data)
	for _, want := range []string{
		"# FDVIB thermochemistry\n",
		"# model: local_harmonic\n",
		"# low_frequency_model: frequency_floor\n",
		"# frequency_floor_cm1: 100\n",
		"# imaginary_modes_excluded: 0\n",
		"# zero_modes_excluded: 0\n",
		"# positive_modes_used: 9\n",
		"# modes_floored: 0\n",
		"# units: T=K energies=eV entropy=eV/K\n",
		"# units: T=K energies=kcal/mol entropy=kcal/mol/K\n",
		"# units: T=K energies=kJ/mol entropy=kJ/mol/K\n",
		"  " + "298.150",
	} {
		if !strings.Contains(text, want) {
			t.Errorf("missing %q in thermo.dat:\n%s", want, text)
		}
	}
}

func TestThermoValidation(t *testing.T) {
	dir := buildResults(t)
	thermoIn := filepath.Join(filepath.Dir(filepath.Dir(dir)), "thermo.in")
	cases := []string{
		"model = bogus\ntemperature_k = 298.15\n",
		"model = local_harmonic\npressure_atm = 1.0\ntemperature_k = 298.15\n",
		"model = gas_rrho\ntemperature_k = 298.15\nzero_tolerance_cm1 = 1.0\n",
		"model = local_harmonic\ntemperature_k = 0\n",
		"model = local_harmonic\ntemperature_k = 298.15\nlow_frequency_model = bogus\n",
		"model = gas_rrho\ntemperature_k = 298.15\nlow_frequency_model = frequency_floor\npressure_atm = 1\nsymmetry_number = 1\n",
		"model = local_harmonic\ntemperature_k = 298.15\nzero_tolerance_cm1 = -1\n",
		"model = local_harmonic\ntemperature_k = 298.15\nlow_frequency_model = frequency_floor\nfrequency_floor_cm1 = 0\n",
		"model = local_harmonic\ntemperature_k = 298.15\nmystery_key = 1\n",
	}
	for i, content := range cases {
		if err := os.WriteFile(thermoIn, []byte(content), 0o644); err != nil {
			t.Fatal(err)
		}
		if err := Thermo(dir, thermoIn); err == nil {
			t.Errorf("case %d should fail: %q", i, content)
		}
	}
}

func TestThermoLocalHarmonicRejectsGasParams(t *testing.T) {
	dir := buildResults(t)
	thermoIn := filepath.Join(filepath.Dir(filepath.Dir(dir)), "thermo.in")
	content := "model = local_harmonic\ntemperature_k = 298.15\nrotor_type = auto\n"
	os.WriteFile(thermoIn, []byte(content), 0o644)
	if err := Thermo(dir, thermoIn); err == nil {
		t.Error("local_harmonic with rotor_type should fail")
	}
}
