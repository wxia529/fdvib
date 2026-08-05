// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package settings assembles and validates the fdvib.in Settings,
// mirroring settings() and job_name() in common.cpp.
package settings

import (
	"fmt"
	"path/filepath"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
)

// Settings is the validated fdvib.in configuration, like struct Settings.
type Settings struct {
	ConfigPath string
	Root       string
	ScfInput   string
	Workdir    string

	SystemType    string
	PWCommand     string
	DynmatCommand string
	OutputPrefix  string

	Selected []int
	// SelectedAll marks selected_atoms = all; Selected is then empty.
	SelectedAll  bool
	Displacement float64
	Multiplicity int
	// MultiplicityExplicit marks an explicit multiplicity in fdvib.in.
	MultiplicityExplicit bool
	RunDynmat            bool
}

// From loads and validates fdvib.in, like settings().
func From(configPath, rootOverride string) (*Settings, error) {
	absPath, err := filepath.Abs(configPath)
	if err != nil {
		return nil, fmt.Errorf("cannot resolve %s", configPath)
	}
	s := &Settings{ConfigPath: absPath}
	if rootOverride == "" {
		s.Root = filepath.Dir(absPath)
	} else {
		abs, err := filepath.Abs(rootOverride)
		if err != nil {
			return nil, fmt.Errorf("cannot resolve %s", rootOverride)
		}
		s.Root = abs
	}
	c, err := config.Load(s.ConfigPath)
	if err != nil {
		return nil, err
	}
	allowed := map[string]bool{
		"scf_input": true, "outdir": true, "system_type": true,
		"selected_atoms": true, "displacement_angstrom": true,
		"multiplicity": true, "pw_command": true, "prefix": true,
		"run_dynmat": true, "dynmat_command": true,
	}
	if err := c.RequireOnly(allowed, "fdvib.in"); err != nil {
		return nil, err
	}
	s.ScfInput = filepath.Join(s.Root, c.Get("scf_input", "scf.in"))
	s.Workdir = filepath.Join(s.Root, c.Get("outdir", "fdvib"))
	s.SystemType = strings.ToLower(c.Get("system_type", "local"))
	s.Displacement, err = c.Real("displacement_angstrom", 0.01)
	if err != nil {
		return nil, err
	}
	s.Multiplicity, err = c.Integer("multiplicity", 1)
	if err != nil {
		return nil, err
	}
	s.MultiplicityExplicit = c.Has("multiplicity")
	s.PWCommand = c.Get("pw_command", "pw.x")
	s.DynmatCommand = c.Get("dynmat_command", "dynmat.x")
	s.OutputPrefix = c.Get("prefix", "system")
	s.RunDynmat, err = config.ParseLogical(c, "run_dynmat", true)
	if err != nil {
		return nil, err
	}
	if s.Displacement <= 0 {
		return nil, fmt.Errorf("displacement_angstrom must be positive")
	}
	if s.Multiplicity < 1 {
		return nil, fmt.Errorf("multiplicity must be positive")
	}
	if s.OutputPrefix == "" || filepath.Base(s.OutputPrefix) != s.OutputPrefix ||
		s.OutputPrefix == "." || s.OutputPrefix == ".." {
		return nil, fmt.Errorf("prefix must be a non-empty filename prefix without directories")
	}
	if strings.Trim(s.PWCommand, " \t\r\n") == "" {
		return nil, fmt.Errorf("pw_command must not be empty")
	}
	if s.RunDynmat && strings.Trim(s.DynmatCommand, " \t\r\n") == "" {
		return nil, fmt.Errorf("dynmat_command must not be empty")
	}
	if s.SystemType != "gas" && s.SystemType != "local" {
		return nil, fmt.Errorf("system_type must be gas or local")
	}
	atoms := strings.ToLower(c.Get("selected_atoms", ""))
	s.SelectedAll = atoms == "all"
	if s.SelectedAll {
		s.Selected = nil
	} else {
		s.Selected, err = config.ParseIntList(atoms)
		if err != nil {
			return nil, err
		}
	}
	return s, nil
}

// JobName builds a displacement job id like disp_0001_x_p, mirroring
// job_name(): axis is 0..2 (x/y/z) and sign > 0 maps to 'p', else 'm'.
func JobName(atom1, axis int, sign int) string {
	xyz := []string{"x", "y", "z"}
	letter := "m"
	if sign > 0 {
		letter = "p"
	}
	return fmt.Sprintf("disp_%04d_%s_%s", atom1, xyz[axis], letter)
}
