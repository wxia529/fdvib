// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package export

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
	modeselect "github.com/wxia529/fdvib/internal/modes"
	"github.com/wxia529/fdvib/internal/qevibration"
	"github.com/wxia529/fdvib/internal/results"
	"github.com/wxia529/fdvib/internal/units"
)

const (
	evToKcalMol = 23.06054783061903
	evToKjMol   = 96.48533212331002
)

// vibrationThermo computes harmonic vibrational thermochemistry, like
// vibration_thermo().
func vibrationThermo(modes []qevibration.Mode, T float64, model string, floor, zeroTol float64) (zpe, u, s, f float64, imag, zero, floored, used int) {
	for _, m := range modes {
		nu := m.Freq
		if math.Abs(nu) < zeroTol {
			zero++
			continue
		}
		if nu < 0 {
			imag++
			continue
		}
		if model == "frequency_floor" && nu < floor {
			nu = floor
			floored++
		}
		used++
		e := nu * units.CMToEV
		x := e / (units.KBEV * T)
		zpe += 0.5 * e
		u += 0.5*e + e/math.Expm1(x)
		s += units.KBEV * (x/math.Expm1(x) - math.Log1p(-math.Exp(-x)))
	}
	f = u - T*s
	return
}

// gasVibrationalModes removes the rigid-body modes and rejects non-positive
// vibrational frequencies, like gas_vibrational_modes().
func gasVibrationalModes(modes []qevibration.Mode, rotorType string) (vibrations []qevibration.Mode, largestRigid float64, err error) {
	rigidDof := 3
	if rotorType == "linear" {
		rigidDof = 5
	} else if rotorType == "nonlinear" {
		rigidDof = 6
	}
	if len(modes) < rigidDof {
		return nil, 0, fmt.Errorf("not enough normal modes for gas RRHO rigid-body projection")
	}
	order := make([]int, len(modes))
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(a, b int) bool {
		return math.Abs(modes[order[a]].Freq) < math.Abs(modes[order[b]].Freq)
	})
	rigid := make([]bool, len(modes))
	largestRigid = 0.0
	for i := 0; i < rigidDof; i++ {
		rigid[order[i]] = true
		if math.Abs(modes[order[i]].Freq) > largestRigid {
			largestRigid = math.Abs(modes[order[i]].Freq)
		}
	}
	vibrations = make([]qevibration.Mode, 0, len(modes)-rigidDof)
	for i := range modes {
		if rigid[i] {
			continue
		}
		if modes[i].Freq <= 0.0 {
			return nil, 0, fmt.Errorf("gas RRHO has an imaginary/non-positive vibrational mode at %s cm^-1; optimize the geometry before thermochemistry",
				config.FormatGeneral(modes[i].Freq, 6))
		}
		vibrations = append(vibrations, modes[i])
	}
	return vibrations, largestRigid, nil
}

// Thermo runs the `thermo` command, like thermo() in thermo.cpp.
func Thermo(resultsDir, thermoInput string) error {
	c, err := config.Load(thermoInput)
	if err != nil {
		return err
	}
	if err := c.RequireOnly(map[string]bool{
		"model": true, "temperature_k": true, "pressure_atm": true,
		"symmetry_number": true, "electronic_degeneracy": true, "rotor_type": true,
		"low_frequency_model": true, "frequency_floor_cm1": true,
		"zero_tolerance_cm1": true,
	}, thermoInput); err != nil {
		return err
	}
	model := strings.ToLower(c.Get("model", ""))
	if model != "gas_rrho" && model != "local_harmonic" {
		return fmt.Errorf("thermo model must be gas_rrho or local_harmonic")
	}
	if model == "local_harmonic" &&
		(c.Has("pressure_atm") || c.Has("symmetry_number") ||
			c.Has("electronic_degeneracy") || c.Has("rotor_type")) {
		return fmt.Errorf("local_harmonic must not contain gas-only thermochemistry parameters")
	}
	if model == "gas_rrho" && c.Has("zero_tolerance_cm1") {
		return fmt.Errorf("gas_rrho must not contain zero_tolerance_cm1; rigid-body modes are excluded by molecular degrees of freedom")
	}
	if _, err := os.Stat(filepath.Join(resultsDir, "dynmat.in")); err == nil {
		d, err := config.LoadQENamelist(filepath.Join(resultsDir, "dynmat.in"))
		if err != nil {
			return err
		}
		asr := strings.ToLower(d.Get("asr", "no"))
		remove := strings.ToLower(d.Get("remove_interaction_blocks", ".false."))
		removesBlocks := remove == ".true." || remove == "true" || remove == "t"
		if asr != "no" {
			return fmt.Errorf("thermochemistry requires asr='no' in dynmat.in")
		}
		if model == "local_harmonic" && !removesBlocks {
			return fmt.Errorf("local_harmonic requires remove_interaction_blocks=.true. in dynmat.in")
		}
		if model == "gas_rrho" && removesBlocks {
			return fmt.Errorf("gas_rrho requires remove_interaction_blocks=.false. in dynmat.in")
		}
	}
	T, err := c.Real("temperature_k", -1)
	if err != nil {
		return err
	}
	if T <= 0 {
		return fmt.Errorf("temperature_k must be > 0 (got %g)", T)
	}
	low := strings.ToLower(c.Get("low_frequency_model", "harmonic"))
	if low != "harmonic" && low != "frequency_floor" {
		return fmt.Errorf("invalid low_frequency_model")
	}
	if model == "gas_rrho" && low != "harmonic" {
		return fmt.Errorf("gas_rrho requires low_frequency_model='harmonic'; rigid-body modes are excluded by molecular degrees of freedom")
	}
	floor, err := c.Real("frequency_floor_cm1", 100)
	if err != nil {
		return err
	}
	zt, err := c.Real("zero_tolerance_cm1", 1)
	if err != nil {
		return err
	}
	if zt < 0.0 {
		return fmt.Errorf("zero_tolerance_cm1 must be non-negative (got %g)", zt)
	}
	if low == "frequency_floor" && floor <= 0.0 {
		return fmt.Errorf("frequency_floor_cm1 must be positive (got %g)", floor)
	}
	dyn, freq, err := results.ResultFiles(resultsDir, "Thermochemistry")
	if err != nil {
		return err
	}
	g, err := qevibration.ReadQeDynGeometry(dyn)
	if err != nil {
		return err
	}
	modes, err := qevibration.ReadQeDynmatModes(freq, len(g.Masses))
	if err != nil {
		return err
	}

	rotorType := "nonlinear"
	var I [3]float64
	mtot := 0.0
	if model == "gas_rrho" {
		for _, m := range g.Masses {
			mtot += m
		}
		var com config.Vec3
		for i := range g.Masses {
			for k := 0; k < 3; k++ {
				com[k] += g.Masses[i] * g.RBohr[i][k] / mtot
			}
		}
		var d, off config.Vec3
		for i := range g.Masses {
			m := g.Masses[i] * units.AmuKg
			var x [3]float64
			for k := 0; k < 3; k++ {
				x[k] = (g.RBohr[i][k] - com[k]) * units.BohrToAng * 1e-10
			}
			d[0] += m * (x[1]*x[1] + x[2]*x[2])
			d[1] += m * (x[0]*x[0] + x[2]*x[2])
			d[2] += m * (x[0]*x[0] + x[1]*x[1])
			off[0] -= m * x[0] * x[1]
			off[1] -= m * x[0] * x[2]
			off[2] -= m * x[1] * x[2]
		}
		I = modeselect.SymmetricEigenvalues(d, off)
		rotorType = strings.ToLower(c.Get("rotor_type", "auto"))
		if rotorType == "auto" {
			if len(g.Masses) == 1 {
				rotorType = "atom"
			} else if I[0]/math.Max(I[2], 1e-300) < 1e-6 {
				rotorType = "linear"
			} else {
				rotorType = "nonlinear"
			}
		}
		if rotorType != "atom" && rotorType != "linear" && rotorType != "nonlinear" {
			return fmt.Errorf("rotor_type must be auto, atom, linear, nonlinear")
		}
		if (rotorType == "atom") != (len(g.Masses) == 1) {
			return fmt.Errorf("rotor_type='atom' is valid only for a monatomic species")
		}
		if rotorType == "nonlinear" && len(g.Masses) < 3 {
			return fmt.Errorf("molecule with fewer than three atoms cannot be a nonlinear rotor")
		}
	}

	largestRigid := 0.0
	var thermoModes []qevibration.Mode
	if model == "gas_rrho" {
		thermoModes, largestRigid, err = gasVibrationalModes(modes, rotorType)
		if err != nil {
			return err
		}
	} else {
		thermoModes = modes
	}
	zeroTol := zt
	if model == "gas_rrho" {
		zeroTol = 0.0
	}
	zpe, vibU, vibS, vibF, imag, zero, floored, used :=
		vibrationThermo(thermoModes, T, low, floor, zeroTol)
	hcorr, stotal, gcorr := vibU, vibS, vibF
	strans, srot, selec, htrans, urot := 0.0, 0.0, 0.0, 0.0, 0.0
	if model == "gas_rrho" {
		if !c.Has("pressure_atm") || !c.Has("symmetry_number") {
			return fmt.Errorf("gas_rrho requires explicit pressure_atm and symmetry_number")
		}
		metadata, err := results.ReadMetadata(resultsDir, true)
		if err != nil {
			return err
		}
		if metadata.ModeSelection != "gas" {
			return fmt.Errorf("gas_rrho requires mode_selection=gas in metadata.dat")
		}
		patm, err := c.Real("pressure_atm", -1)
		if err != nil {
			return err
		}
		sigma, err := c.Integer("symmetry_number", 0)
		if err != nil {
			return err
		}
		mult := metadata.Multiplicity
		if patm <= 0 || sigma <= 0 || mult <= 0 {
			return fmt.Errorf("gas pressure, symmetry, and multiplicity must be positive")
		}
		qtrans := math.Pow(2*units.PI*(mtot*units.AmuKg)*units.KBSI*T/(units.HSI*units.HSI), 1.5) *
			(units.KBSI * T / (patm * units.AtmPa))
		strans = units.KBEV * (math.Log(qtrans) + 2.5)
		htrans = 2.5 * units.KBEV * T
		if rotorType == "linear" {
			moment := math.Max(I[1], I[2])
			if !(moment > 0.0) {
				return fmt.Errorf("invalid linear-molecule moment of inertia")
			}
			qr := 8 * units.PI * units.PI * moment * units.KBSI * T / (float64(sigma) * units.HSI * units.HSI)
			srot = units.KBEV * (math.Log(qr) + 1)
			urot = units.KBEV * T
		} else if rotorType == "nonlinear" {
			if !(I[0] > 0.0 && I[1] > 0.0 && I[2] > 0.0) {
				return fmt.Errorf("invalid nonlinear-molecule moments of inertia")
			}
			qr := math.Sqrt(units.PI) / float64(sigma) *
				math.Pow(8*units.PI*units.PI*units.KBSI*T/(units.HSI*units.HSI), 1.5) *
				math.Sqrt(I[0]*I[1]*I[2])
			srot = units.KBEV * (math.Log(qr) + 1.5)
			urot = 1.5 * units.KBEV * T
		}
		degeneracy := mult
		ed := strings.ToLower(c.Get("electronic_degeneracy", "auto"))
		if ed != "auto" {
			value, err := config.ParseNumber(ed)
			if err != nil {
				return err
			}
			rounded := math.Round(value)
			if math.Abs(value-rounded) > 1.0e-10 {
				return fmt.Errorf("electronic_degeneracy must be an integer or auto")
			}
			degeneracy = int(rounded)
		}
		if degeneracy <= 0 {
			return fmt.Errorf("electronic_degeneracy must be positive or auto (got %d)", degeneracy)
		}
		selec = units.KBEV * math.Log(float64(degeneracy))
		hcorr = vibU + htrans + urot
		stotal = vibS + strans + srot + selec
		gcorr = hcorr - T*stotal
	}

	var b strings.Builder
	fmt.Fprintf(&b, "# FDVIB thermochemistry\n# model: %s\n# low_frequency_model: %s\n", model, low)
	if low == "frequency_floor" {
		fmt.Fprintf(&b, "# frequency_floor_cm1: %s\n", config.FormatGeneral(floor, 6))
	}
	if model == "gas_rrho" {
		rigidDof := 3
		if rotorType == "linear" {
			rigidDof = 5
		} else if rotorType == "nonlinear" {
			rigidDof = 6
		}
		fmt.Fprintf(&b, "# rotor_type: %s\n# rigid_body_modes_excluded: %d\n# max_rigid_body_frequency_cm1: %s\n# expected_vibrational_modes: %d\n",
			rotorType, rigidDof, config.FormatGeneral(largestRigid, 6), len(thermoModes))
	}
	fmt.Fprintf(&b, "# imaginary_modes_excluded: %d\n# zero_modes_excluded: %d\n# positive_modes_used: %d\n# modes_floored: %d\n",
		imag, zero, used, floored)
	b.WriteString("# units: T=K energies=eV entropy=eV/K\n")
	if model == "local_harmonic" {
		b.WriteString("# " + fmt.Sprintf("%11s", "T/K") + fmt.Sprintf("%16s", "ZPE/eV") + fmt.Sprintf("%16s", "U_vib/eV") +
			fmt.Sprintf("%19s", "S_vib/eV_K") + fmt.Sprintf("%16s", "TS_vib/eV") + fmt.Sprintf("%16s", "F_vib/eV") + "\n")
		fmt.Fprintf(&b, "  %s%s%s%s%s%s\n",
			config.FixedField(T, 3, 11),
			config.FixedField(zpe, 10, 16),
			config.FixedField(vibU, 10, 16),
			config.FixedField(vibS, 12, 19),
			config.FixedField(T*vibS, 10, 16),
			config.FixedField(vibF, 10, 16))
	} else {
		b.WriteString("# " + fmt.Sprintf("%11s", "T/K") + fmt.Sprintf("%15s", "ZPE/eV") + fmt.Sprintf("%15s", "U_vib/eV") +
			fmt.Sprintf("%15s", "H_trans/eV") + fmt.Sprintf("%15s", "U_rot/eV") +
			fmt.Sprintf("%18s", "S_trans/eV_K") + fmt.Sprintf("%18s", "S_rot/eV_K") +
			fmt.Sprintf("%18s", "S_vib/eV_K") + fmt.Sprintf("%18s", "S_elec/eV_K") +
			fmt.Sprintf("%15s", "H_corr/eV") + fmt.Sprintf("%15s", "G_corr/eV") + "\n")
		fmt.Fprintf(&b, "  %s%s%s%s%s%s%s%s%s%s%s\n",
			config.FixedField(T, 3, 11),
			config.FixedField(zpe, 10, 15),
			config.FixedField(vibU, 10, 15),
			config.FixedField(htrans, 10, 15),
			config.FixedField(urot, 10, 15),
			config.FixedField(strans, 12, 18),
			config.FixedField(srot, 12, 18),
			config.FixedField(vibS, 12, 18),
			config.FixedField(selec, 12, 18),
			config.FixedField(hcorr, 10, 15),
			config.FixedField(gcorr, 10, 15))
	}
	appendMolarTable := func(scale float64, energyUnit string) {
		b.WriteString("\n# units: T=K energies=" + energyUnit + " entropy=" + energyUnit + "/K\n")
		if model == "local_harmonic" {
			b.WriteString("# " + fmt.Sprintf("%11s", "T/K") + fmt.Sprintf("%16s", "ZPE") + fmt.Sprintf("%16s", "U_vib") +
				fmt.Sprintf("%19s", "S_vib") + fmt.Sprintf("%16s", "TS_vib") + fmt.Sprintf("%16s", "F_vib") + "\n")
			fmt.Fprintf(&b, "  %s%s%s%s%s%s\n",
				config.FixedField(T, 3, 11),
				config.FixedField(zpe*scale, 10, 16),
				config.FixedField(vibU*scale, 10, 16),
				config.FixedField(vibS*scale, 12, 19),
				config.FixedField(T*vibS*scale, 10, 16),
				config.FixedField(vibF*scale, 10, 16))
		} else {
			b.WriteString("# " + fmt.Sprintf("%11s", "T/K") + fmt.Sprintf("%15s", "ZPE") + fmt.Sprintf("%15s", "U_vib") +
				fmt.Sprintf("%15s", "H_trans") + fmt.Sprintf("%15s", "U_rot") +
				fmt.Sprintf("%18s", "S_trans") + fmt.Sprintf("%18s", "S_rot") +
				fmt.Sprintf("%18s", "S_vib") + fmt.Sprintf("%18s", "S_elec") +
				fmt.Sprintf("%15s", "H_corr") + fmt.Sprintf("%15s", "G_corr") + "\n")
			fmt.Fprintf(&b, "  %s%s%s%s%s%s%s%s%s%s%s\n",
				config.FixedField(T, 3, 11),
				config.FixedField(zpe*scale, 10, 15),
				config.FixedField(vibU*scale, 10, 15),
				config.FixedField(htrans*scale, 10, 15),
				config.FixedField(urot*scale, 10, 15),
				config.FixedField(strans*scale, 12, 18),
				config.FixedField(srot*scale, 12, 18),
				config.FixedField(vibS*scale, 12, 18),
				config.FixedField(selec*scale, 12, 18),
				config.FixedField(hcorr*scale, 10, 15),
				config.FixedField(gcorr*scale, 10, 15))
		}
	}
	appendMolarTable(evToKcalMol, "kcal/mol")
	appendMolarTable(evToKjMol, "kJ/mol")
	out := b.String()
	if err := config.WriteText(filepath.Join(resultsDir, "thermo.dat"), out); err != nil {
		return err
	}
	fmt.Print(out)
	fmt.Printf("Written %s\n", config.DisplayPath(filepath.Join(resultsDir, "thermo.dat")))
	return nil
}
