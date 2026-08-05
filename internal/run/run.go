// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package run orchestrates the five-stage calculation (reference SCF,
// displacements, Hessian analysis, dynmat.x) with locking, seeding,
// recovery, and integrity digests, mirroring prepare_run.cpp.
package run

import (
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"syscall"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/hessian"
	"github.com/wxia529/fdvib/internal/process"
	"github.com/wxia529/fdvib/internal/qeinput"
	"github.com/wxia529/fdvib/internal/qeoutput"
	"github.com/wxia529/fdvib/internal/qevibration"
	"github.com/wxia529/fdvib/internal/settings"
	"github.com/wxia529/fdvib/internal/state"
)

// ReferenceSeed carries the converged reference charge density (and optional
// PAW data) plus their digests, like struct ReferenceSeed.
type ReferenceSeed struct {
	Density       string
	Paw           string
	DensityDigest string
	PawDigest     string
}

var (
	nspinRe         = regexp.MustCompile(`(?i)\bnspin\s*=\s*(\d+)`)
	magnetizationRe = regexp.MustCompile(`(?i)\btot_magnetization\s*=\s*([-+0-9.EeDd]+)`)
)

// Calculate runs the full calculation, like calculate() in prepare_run.cpp.
func Calculate(s *settings.Settings) error {
	unlock, err := acquireLock(s.Workdir)
	if err != nil {
		return err
	}
	defer unlock()

	q, err := qeinput.ParseQeInput(s.ScfInput)
	if err != nil {
		return err
	}
	selected, err := selectedAtoms(s, q.Nat)
	if err != nil {
		return err
	}
	if s.SystemType == "gas" {
		nspin := 1
		if m := nspinRe.FindStringSubmatch(q.CleanText); m != nil {
			nspin, _ = strconv.Atoi(m[1])
		}
		magnetization := 0.0
		hasMagnetization := false
		if m := magnetizationRe.FindStringSubmatch(q.CleanText); m != nil {
			x, err := config.Number(m[1])
			if err != nil {
				return err
			}
			magnetization = x
			hasMagnetization = true
		}
		if s.Multiplicity == 1 {
			if nspin != 1 && !(nspin == 2 && hasMagnetization && math.Abs(magnetization) < 1e-8) {
				return fmt.Errorf("Gas singlet requires nspin=1, or nspin=2 with tot_magnetization=0")
			}
		} else if nspin != 2 || !hasMagnetization ||
			math.Abs(magnetization-float64(s.Multiplicity-1)) > 1e-8 {
			return fmt.Errorf("Gas multiplicity requires nspin=2 and tot_magnetization=multiplicity-1")
		}
	}
	if err := initializeDataset(s, q, selected); err != nil {
		return err
	}
	reference, err := ensureReference(s, q)
	if err != nil {
		return err
	}
	if err := runDisplacements(s, q, selected, reference); err != nil {
		return err
	}
	if err := ensureAnalysis(s); err != nil {
		return err
	}
	if err := ensureDynmat(s); err != nil {
		return err
	}
	fmt.Println("FDVIB calculation completed")
	return nil
}

// acquireLock takes an exclusive flock on <workdir>.lock and returns the
// release function, like CalculationLock.
func acquireLock(workdir string) (func(), error) {
	lockPath := workdir + ".lock"
	f, err := os.OpenFile(lockPath, os.O_CREATE|os.O_RDWR, 0o644)
	if err != nil {
		return nil, fmt.Errorf("Cannot open calculation lock: %s", lockPath)
	}
	if err := syscall.Flock(int(f.Fd()), syscall.LOCK_EX|syscall.LOCK_NB); err != nil {
		f.Close()
		return nil, fmt.Errorf("Calculation is already running: %s", workdir)
	}
	if err := f.Truncate(0); err == nil {
		fmt.Fprintf(f, "pid=%d\n", os.Getpid())
	}
	return func() {
		_ = syscall.Flock(int(f.Fd()), syscall.LOCK_UN)
		f.Close()
	}, nil
}

func metadataText(s *settings.Settings, output string) (string, error) {
	energy, err := qeoutput.ReadTotalEnergyHartree(output)
	if err != nil {
		return "", err
	}
	var b strings.Builder
	fmt.Fprintf(&b, "program = qe\n")
	fmt.Fprintf(&b, "electronic_energy_hartree = %s\n", config.FormatSci(energy, 16))
	fmt.Fprintf(&b, "multiplicity = %d\n", s.Multiplicity)
	fmt.Fprintf(&b, "mode_selection = %s\n", s.SystemType)
	b.WriteString("selected_atoms = ")
	if s.SelectedAll {
		b.WriteString("all")
	} else {
		for i, atom := range s.Selected {
			if i > 0 {
				b.WriteString(",")
			}
			fmt.Fprintf(&b, "%d", atom)
		}
	}
	b.WriteString("\n")
	return b.String(), nil
}

// densityFile finds the single non-empty charge density of an attempt,
// like density_file().
func densityFile(attempt, prefix string) (string, error) {
	save := filepath.Join(attempt, "out", prefix+".save")
	var found []string
	for _, name := range []string{"charge-density.dat", "charge-density.hdf5"} {
		p := filepath.Join(save, name)
		if fi, err := os.Stat(p); err == nil && fi.Mode().IsRegular() && fi.Size() > 0 {
			found = append(found, p)
		}
	}
	if len(found) != 1 {
		return "", fmt.Errorf("Expected exactly one non-empty charge-density.dat or charge-density.hdf5 in %s", save)
	}
	return found[0], nil
}

// removeDisplacementDensity deletes the copied reference density from a
// completed displacement attempt, like remove_displacement_density().
func removeDisplacementDensity(attempt, prefix string, referenceDensity string) error {
	copied := filepath.Join(attempt, "out", prefix+".save", filepath.Base(referenceDensity))
	err := os.Remove(copied)
	if err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("Cannot remove completed displacement charge density %s: %v", copied, err)
	}
	if _, statErr := os.Stat(copied); statErr == nil {
		return fmt.Errorf("Cannot remove completed displacement charge density: %s", copied)
	}
	return nil
}

// selectedAtoms validates and resolves the selected atom list, like
// selected_atoms() in prepare_run.cpp.
func selectedAtoms(s *settings.Settings, nat int) ([]int, error) {
	selected := s.Selected
	if s.SystemType == "gas" {
		if !s.MultiplicityExplicit || s.Multiplicity < 1 {
			return nil, fmt.Errorf("gas requires an explicit positive multiplicity in fdvib.in")
		}
		if !s.SelectedAll {
			return nil, fmt.Errorf("gas requires selected_atoms='all'")
		}
		selected = make([]int, nat)
		for i := range selected {
			selected[i] = i + 1
		}
	}
	if len(selected) == 0 {
		return nil, fmt.Errorf("selected_atoms is empty")
	}
	unique := make(map[int]bool)
	for _, atom := range selected {
		if atom < 1 || atom > nat || unique[atom] {
			return nil, fmt.Errorf("Invalid/duplicate selected atom")
		}
		unique[atom] = true
	}
	return selected, nil
}

func datasetDescription(s *settings.Settings, q *qeinput.QEInput, selected []int) (string, error) {
	scfDigest, err := state.FileDigest(s.ScfInput)
	if err != nil {
		return "", err
	}
	var b strings.Builder
	b.WriteString("format=2\n")
	fmt.Fprintf(&b, "scf_digest=%s\n", scfDigest)
	fmt.Fprintf(&b, "system_type=%s\n", s.SystemType)
	fmt.Fprintf(&b, "displacement_angstrom=%s\n", config.FormatGeneral(s.Displacement, 17))
	fmt.Fprintf(&b, "multiplicity=%d\n", s.Multiplicity)
	fmt.Fprintf(&b, "prefix=%s\n", s.OutputPrefix)
	fmt.Fprintf(&b, "qe_prefix=%s\nselected_atoms=", q.Prefix)
	for _, atom := range selected {
		fmt.Fprintf(&b, "%d,", atom)
	}
	b.WriteString("\n")
	return b.String(), nil
}

func initializeDataset(s *settings.Settings, q *qeinput.QEInput, selected []int) error {
	stateDir := filepath.Join(s.Workdir, "state")
	wanted, err := datasetDescription(s, q, selected)
	if err != nil {
		return err
	}
	if fi, err := os.Stat(s.Workdir); err == nil && fi.IsDir() {
		if _, err := os.Stat(filepath.Join(stateDir, "dataset.state")); err != nil {
			entries, readErr := os.ReadDir(s.Workdir)
			if readErr == nil && len(entries) > 0 {
				return fmt.Errorf("Refusing to use non-empty outdir without FDVIB state metadata: %s", s.Workdir)
			}
		}
	}
	if err := os.MkdirAll(stateDir, 0o755); err != nil {
		return fmt.Errorf("Cannot create %s", stateDir)
	}
	dataset := filepath.Join(stateDir, "dataset.state")
	if _, err := os.Stat(dataset); err == nil {
		existing, readErr := config.ReadText(dataset)
		if readErr != nil {
			return readErr
		}
		if existing != wanted {
			return fmt.Errorf("Dataset differs from the existing calculation; use a different outdir")
		}
	} else if err := config.WriteText(dataset, wanted); err != nil {
		return err
	}
	scfReference := filepath.Join(s.Workdir, "scf.in.reference")
	if _, err := os.Stat(scfReference); err != nil {
		if err := copyFile(s.ScfInput, scfReference); err != nil {
			return err
		}
	} else {
		ref, err := config.ReadText(scfReference)
		if err != nil {
			return err
		}
		input, err := config.ReadText(s.ScfInput)
		if err != nil {
			return err
		}
		if ref != input {
			return fmt.Errorf("scf.in differs from the existing calculation")
		}
	}
	var cfg strings.Builder
	cfg.WriteString("scf_input = scf.in\noutdir = fdvib\n")
	fmt.Fprintf(&cfg, "system_type = %s\nselected_atoms = ", s.SystemType)
	if s.SelectedAll {
		cfg.WriteString("all")
	} else {
		for i, atom := range selected {
			if i > 0 {
				cfg.WriteString(",")
			}
			fmt.Fprintf(&cfg, "%d", atom)
		}
	}
	cfg.WriteString("\n")
	fmt.Fprintf(&cfg, "displacement_angstrom = %s\n", config.FormatGeneral(s.Displacement, 17))
	fmt.Fprintf(&cfg, "multiplicity = %d\nprefix = %s\n", s.Multiplicity, s.OutputPrefix)
	cfg.WriteString("pw_command = pw.x\nrun_dynmat = false\ndynmat_command = dynmat.x\n")
	configReference := filepath.Join(s.Workdir, "fdvib.in.reference")
	if _, err := os.Stat(configReference); err != nil {
		if err := config.WriteText(configReference, cfg.String()); err != nil {
			return err
		}
	} else {
		existing, err := config.ReadText(configReference)
		if err != nil {
			return err
		}
		if existing != cfg.String() {
			return fmt.Errorf("FDVIB dataset snapshot is missing or modified: %s", configReference)
		}
	}
	return nil
}

func copyFile(src, dst string) error {
	in, err := os.Open(src)
	if err != nil {
		return fmt.Errorf("Cannot read %s", src)
	}
	defer in.Close()
	fi, err := in.Stat()
	if err != nil {
		return fmt.Errorf("Cannot read %s", src)
	}
	if err := os.MkdirAll(filepath.Dir(dst), 0o755); err != nil {
		return fmt.Errorf("Cannot write %s", dst)
	}
	// fs::copy_file preserves the source permission bits.
	out, err := os.OpenFile(dst, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, fi.Mode().Perm())
	if err != nil {
		return fmt.Errorf("Cannot write %s", dst)
	}
	defer out.Close()
	if _, err := io.Copy(out, in); err != nil {
		return fmt.Errorf("Cannot write %s", dst)
	}
	return nil
}

func ensureReference(s *settings.Settings, q *qeinput.QEInput) (*ReferenceSeed, error) {
	calculations := filepath.Join(s.Workdir, "calculations")
	marker := filepath.Join(s.Workdir, "state", "init_scf.complete")
	if _, err := os.Stat(marker); err == nil {
		text, err := config.ReadText(marker)
		if err != nil {
			return nil, err
		}
		fields := strings.Fields(text)
		if len(fields) != 5 || !state.IsNumberedName(fields[0], "init_scf") {
			return nil, fmt.Errorf("Invalid reference completion snapshot: %s", marker)
		}
		attemptName, densityName, densityDigest, pawName, pawDigest := fields[0], fields[1], fields[2], fields[3], fields[4]
		attempt := filepath.Join(calculations, attemptName)
		density := filepath.Join(attempt, "out", q.Prefix+".save", densityName)
		paw := filepath.Join(filepath.Dir(density), "paw.txt")
		if err := qeoutput.ValidateQeOutput(filepath.Join(attempt, "scf.out")); err != nil {
			return nil, err
		}
		if !isRegularFile(density) {
			return nil, fmt.Errorf("Completed reference charge density is missing or modified: %s", density)
		}
		d, err := state.FileDigest(density)
		if err != nil {
			return nil, err
		}
		if d != densityDigest {
			return nil, fmt.Errorf("Completed reference charge density is missing or modified: %s", density)
		}
		if pawName == "-" {
			if pawDigest != "-" || fileExists(paw) {
				return nil, fmt.Errorf("Completed reference PAW snapshot is inconsistent: %s", paw)
			}
		} else if pawName != "paw.txt" || !isRegularNonempty(paw) {
			return nil, fmt.Errorf("Completed reference PAW data is missing or modified: %s", paw)
		} else if pd, err := state.FileDigest(paw); err != nil || pd != pawDigest {
			return nil, fmt.Errorf("Completed reference PAW data is missing or modified: %s", paw)
		}
		wantedMetadata, err := metadataText(s, filepath.Join(attempt, "scf.out"))
		if err != nil {
			return nil, err
		}
		metadataPath := filepath.Join(s.Workdir, "metadata.dat")
		if !isRegularFile(metadataPath) {
			return nil, fmt.Errorf("Reference metadata is missing or modified: %s", metadataPath)
		}
		existing, err := config.ReadText(metadataPath)
		if err != nil {
			return nil, err
		}
		if existing != wantedMetadata {
			return nil, fmt.Errorf("Reference metadata is missing or modified: %s", metadataPath)
		}
		fmt.Println("Preserved completed reference SCF")
		seed := &ReferenceSeed{Density: density, DensityDigest: densityDigest, PawDigest: pawDigest}
		if pawName == "paw.txt" {
			seed.Paw = paw
		}
		return seed, nil
	}

	commit := func(calculation string, recovered bool) (*ReferenceSeed, error) {
		output := filepath.Join(calculation, "scf.out")
		if err := qeoutput.ValidateQeOutput(output); err != nil {
			return nil, err
		}
		if _, err := qeoutput.ParseForces(output, q.Nat); err != nil {
			return nil, err
		}
		density, err := densityFile(calculation, q.Prefix)
		if err != nil {
			return nil, err
		}
		paw := filepath.Join(filepath.Dir(density), "paw.txt")
		if fi, err := os.Stat(paw); err == nil && (!fi.Mode().IsRegular() || fi.Size() == 0) {
			return nil, fmt.Errorf("Reference PAW data is not a non-empty regular file: %s", paw)
		}
		wantedMetadata, err := metadataText(s, output)
		if err != nil {
			return nil, err
		}
		if err := config.WriteText(filepath.Join(s.Workdir, "metadata.dat"), wantedMetadata); err != nil {
			return nil, err
		}
		hasPaw := isRegularFile(paw)
		densityDigest, err := state.FileDigest(density)
		if err != nil {
			return nil, err
		}
		pawDigest := ""
		if hasPaw {
			pawDigest, err = state.FileDigest(paw)
			if err != nil {
				return nil, err
			}
		}
		markerText := filepath.Base(calculation) + " " + filepath.Base(density) + " " + densityDigest + " "
		if hasPaw {
			markerText += "paw.txt " + pawDigest + "\n"
		} else {
			markerText += "- -\n"
		}
		if err := config.WriteText(marker, markerText); err != nil {
			return nil, err
		}
		if recovered {
			fmt.Printf("Recovered completed initial SCF from %s\n", filepath.Base(calculation))
		}
		seed := &ReferenceSeed{Density: density, DensityDigest: densityDigest, PawDigest: pawDigest}
		if hasPaw {
			seed.Paw = paw
		}
		return seed, nil
	}
	for _, calculation := range state.NumberedDirectories(calculations, "init_scf") {
		seed, err := commit(calculation, true)
		if err == nil {
			return seed, nil
		}
		// Retain incomplete calculations and try a fresh numbered directory below.
	}
	attempt, err := state.NewNumberedDirectory(calculations, "init_scf")
	if err != nil {
		return nil, err
	}
	input, err := qeinput.ReferenceInput(q, "./out", s.Root, attempt)
	if err != nil {
		return nil, err
	}
	if err := config.WriteText(filepath.Join(attempt, "scf.in"), input); err != nil {
		return nil, err
	}
	cmd := s.PWCommand + " -inp scf.in"
	fmt.Println("Running unperturbed reference SCF")
	rc, err := process.ShellRun(cmd, attempt, filepath.Join(attempt, "scf.out"))
	if err != nil {
		return nil, err
	}
	if rc != 0 {
		return nil, fmt.Errorf("Reference SCF failed with exit code %d", rc)
	}
	return commit(attempt, false)
}

func runDisplacements(s *settings.Settings, q *qeinput.QEInput, selected []int,
	reference *ReferenceSeed) error {
	completed, preserved := 0, 0
	calculations := filepath.Join(s.Workdir, "calculations")
	for _, atom1 := range selected {
		for axis := 0; axis < 3; axis++ {
			for _, sign := range []int{1, -1} {
				id := settings.JobName(atom1, axis, sign)
				marker := filepath.Join(s.Workdir, "state", id+".complete")
				if _, err := os.Stat(marker); err == nil {
					text, err := config.ReadText(marker)
					if err != nil {
						return err
					}
					fields := strings.Fields(text)
					if len(fields) != 2 || !state.IsNumberedName(fields[0], id) {
						return fmt.Errorf("Invalid displacement completion snapshot: %s", marker)
					}
					savedAttempt, digest := fields[0], fields[1]
					attempt := filepath.Join(calculations, savedAttempt)
					output := filepath.Join(attempt, "pw.out")
					forces := filepath.Join(attempt, "forces.dat")
					if _, err := qeoutput.ParseForces(output, q.Nat); err != nil {
						return err
					}
					if !isRegularFile(forces) {
						return fmt.Errorf("Completed force data is missing or modified: %s", forces)
					}
					d, err := state.FileDigest(forces)
					if err != nil {
						return err
					}
					if d != digest {
						return fmt.Errorf("Completed force data is missing or modified: %s", forces)
					}
					if err := removeDisplacementDensity(attempt, q.Prefix, reference.Density); err != nil {
						return err
					}
					preserved++
					continue
				}
				commit := func(calculation string) error {
					output := filepath.Join(calculation, "pw.out")
					forces := filepath.Join(calculation, "forces.dat")
					parsed, err := qeoutput.ParseForces(output, q.Nat)
					if err != nil {
						return err
					}
					if err := qeoutput.WriteForces(forces, parsed, output); err != nil {
						return err
					}
					digest, err := state.FileDigest(forces)
					if err != nil {
						return err
					}
					return config.WriteText(marker, filepath.Base(calculation)+" "+digest+"\n")
				}
				var recovered string
				for _, calculation := range state.NumberedDirectories(calculations, id) {
					if err := commit(calculation); err == nil {
						recovered = calculation
						break
					}
					// Retain incomplete calculations and try a fresh numbered directory below.
				}
				if recovered != "" {
					if err := removeDisplacementDensity(recovered, q.Prefix, reference.Density); err != nil {
						return err
					}
					fmt.Printf("Recovered completed %s from %s\n", id, filepath.Base(recovered))
					completed++
					continue
				}
				attempt, err := state.NewNumberedDirectory(calculations, id)
				if err != nil {
					return err
				}
				input, err := qeinput.DisplacedInput(q, atom1-1, axis, float64(sign)*s.Displacement,
					"./out", s.Root, attempt)
				if err != nil {
					return err
				}
				if err := config.WriteText(filepath.Join(attempt, "pw.in"), input); err != nil {
					return err
				}
				seeded := filepath.Join(attempt, "out", q.Prefix+".save", filepath.Base(reference.Density))
				if err := os.MkdirAll(filepath.Dir(seeded), 0o755); err != nil {
					return fmt.Errorf("Cannot create %s", filepath.Dir(seeded))
				}
				if err := copyFile(reference.Density, seeded); err != nil {
					return err
				}
				d, err := state.FileDigest(seeded)
				if err != nil {
					return err
				}
				if d != reference.DensityDigest {
					return fmt.Errorf("Copied reference charge density failed verification: %s", seeded)
				}
				if reference.Paw != "" {
					seededPaw := filepath.Join(filepath.Dir(seeded), "paw.txt")
					if err := copyFile(reference.Paw, seededPaw); err != nil {
						return err
					}
					pd, err := state.FileDigest(seededPaw)
					if err != nil {
						return err
					}
					if pd != reference.PawDigest {
						return fmt.Errorf("Copied reference PAW data failed verification: %s", seededPaw)
					}
				}
				cmd := s.PWCommand + " -inp pw.in"
				fmt.Printf("Running %s\n", id)
				rc, err := process.ShellRun(cmd, attempt, filepath.Join(attempt, "pw.out"))
				if err != nil {
					return err
				}
				if rc != 0 {
					return fmt.Errorf("%s failed with exit code %d", id, rc)
				}
				if err := commit(attempt); err != nil {
					return err
				}
				if err := removeDisplacementDensity(attempt, q.Prefix, reference.Density); err != nil {
					return err
				}
				completed++
			}
		}
	}
	fmt.Printf("Completed %d, preserved %d displacement jobs\n", completed, preserved)
	return nil
}

func ensureAnalysis(s *settings.Settings) error {
	marker := filepath.Join(s.Workdir, "state", "analyze.complete")
	resultsDir := filepath.Join(s.Workdir, "results")
	dyn := filepath.Join(resultsDir, s.OutputPrefix+".dynG")
	resultMetadata := filepath.Join(resultsDir, "metadata.dat")
	dynmatIn := filepath.Join(resultsDir, "dynmat.in")
	if _, err := os.Stat(marker); err == nil {
		if !isRegularFile(dyn) || !isRegularFile(dynmatIn) {
			return fmt.Errorf("Completed Hessian results are missing")
		}
		text, err := config.ReadText(marker)
		if err != nil {
			return err
		}
		fields := strings.Fields(text)
		if len(fields) != 3 {
			return fmt.Errorf("Invalid Hessian completion snapshot: %s", marker)
		}
		dynDigest, inputDigest, metadataDigest := fields[0], fields[1], fields[2]
		if !isRegularFile(resultMetadata) {
			return fmt.Errorf("Completed Hessian results were modified")
		}
		ok, err := digestsMatch(dyn, dynDigest)
		if err != nil {
			return err
		}
		ok2, err := digestsMatch(dynmatIn, inputDigest)
		if err != nil {
			return err
		}
		ok3, err := digestsMatch(resultMetadata, metadataDigest)
		if err != nil {
			return err
		}
		if !ok || !ok2 || !ok3 {
			return fmt.Errorf("Completed Hessian results were modified")
		}
		fmt.Println("Preserved completed Hessian analysis")
		return nil
	}
	if _, err := os.Stat(resultsDir); err == nil {
		recovered := false
		if isRegularFile(dyn) && isRegularFile(dynmatIn) && isRegularFile(resultMetadata) {
			if geometry, err := qevibration.ReadQeDynGeometry(dyn); err == nil {
				_ = geometry
				if input, err := config.LoadQENamelist(dynmatIn); err == nil &&
					input.Get("fildyn", "") == filepath.Base(dyn) {
					if meta, err := config.ReadText(resultMetadata); err == nil {
						if workMeta, err := config.ReadText(filepath.Join(s.Workdir, "metadata.dat")); err == nil &&
							meta == workMeta {
							d1, err1 := state.FileDigest(dyn)
							d2, err2 := state.FileDigest(dynmatIn)
							d3, err3 := state.FileDigest(resultMetadata)
							if err1 == nil && err2 == nil && err3 == nil &&
								config.WriteText(marker, d1+" "+d2+" "+d3+"\n") == nil {
								fmt.Println("Recovered completed Hessian analysis")
								recovered = true
							}
						}
					}
				}
			}
		}
		if recovered {
			return nil
		}
		failedRoot := filepath.Join(s.Workdir, "failed")
		failed, err := state.NewNumberedDirectory(failedRoot, "analysis")
		if err != nil {
			return err
		}
		if err := syscall.Rename(resultsDir, failed); err != nil {
			return fmt.Errorf("Cannot move %s to %s", resultsDir, failed)
		}
		fmt.Printf("Preserved incomplete Hessian results in %s\n", config.DisplayPath(failed))
	}
	if err := hessian.Analyze(s); err != nil {
		return err
	}
	d1, err := state.FileDigest(dyn)
	if err != nil {
		return err
	}
	d2, err := state.FileDigest(dynmatIn)
	if err != nil {
		return err
	}
	d3, err := state.FileDigest(resultMetadata)
	if err != nil {
		return err
	}
	return config.WriteText(marker, d1+" "+d2+" "+d3+"\n")
}

func ensureDynmat(s *settings.Settings) error {
	if !s.RunDynmat {
		fmt.Println("dynmat.x was not requested (run_dynmat=.false.)")
		return nil
	}
	stateMarker := filepath.Join(s.Workdir, "state", "dynmat.complete")
	resultsDir := filepath.Join(s.Workdir, "results")
	finalOutput := filepath.Join(resultsDir, "dynmat.out")
	finalFreq := filepath.Join(resultsDir, s.OutputPrefix+".freq.out")
	dyn := filepath.Join(resultsDir, s.OutputPrefix+".dynG")
	calculations := filepath.Join(s.Workdir, "calculations")

	validateModes := func() error {
		geometry, err := qevibration.ReadQeDynGeometry(dyn)
		if err != nil {
			return err
		}
		_, err = qevibration.ReadQeDynmatModes(finalFreq, len(geometry.Masses))
		return err
	}
	validateCalculation := func(calculation string) error {
		calculationDyn := filepath.Join(calculation, filepath.Base(dyn))
		calculationInput := filepath.Join(calculation, "dynmat.in")
		calculationOutput := filepath.Join(calculation, "dynmat.out")
		calculationFreq := filepath.Join(calculation, filepath.Base(finalFreq))
		if !isRegularFile(calculationDyn) || !isRegularFile(calculationInput) {
			return fmt.Errorf("dynmat calculation inputs differ from Hessian results: %s", calculation)
		}
		d1, err := state.FileDigest(calculationDyn)
		if err != nil {
			return err
		}
		d2, err := state.FileDigest(calculationInput)
		if err != nil {
			return err
		}
		d3, err := state.FileDigest(dyn)
		if err != nil {
			return err
		}
		d4, err := state.FileDigest(filepath.Join(resultsDir, "dynmat.in"))
		if err != nil {
			return err
		}
		if d1 != d3 || d2 != d4 {
			return fmt.Errorf("dynmat calculation inputs differ from Hessian results: %s", calculation)
		}
		if err := qeoutput.ValidateQeOutput(calculationOutput); err != nil {
			return err
		}
		if !isRegularNonempty(calculationFreq) {
			return fmt.Errorf("dynmat calculation frequency output is missing: %s", calculationFreq)
		}
		geometry, err := qevibration.ReadQeDynGeometry(calculationDyn)
		if err != nil {
			return err
		}
		_, err = qevibration.ReadQeDynmatModes(calculationFreq, len(geometry.Masses))
		return err
	}

	if _, err := os.Stat(stateMarker); err == nil {
		text, err := config.ReadText(stateMarker)
		if err != nil {
			return err
		}
		fields := strings.Fields(text)
		if len(fields) != 3 || !state.IsNumberedName(fields[0], "dynmat") {
			return fmt.Errorf("Invalid dynmat completion snapshot: %s", stateMarker)
		}
		calculationName, outputDigest, freqDigest := fields[0], fields[1], fields[2]
		calculation := filepath.Join(calculations, calculationName)
		if err := validateCalculation(calculation); err != nil {
			return err
		}
		if err := qeoutput.ValidateQeOutput(finalOutput); err != nil {
			return err
		}
		if !isRegularFile(finalFreq) {
			return fmt.Errorf("Completed dynmat frequency output is missing")
		}
		ok1, err := digestsMatch(finalOutput, outputDigest)
		if err != nil {
			return err
		}
		ok2, err := digestsMatch(finalFreq, freqDigest)
		if err != nil {
			return err
		}
		ok3, err := digestsMatch(filepath.Join(calculation, "dynmat.out"), outputDigest)
		if err != nil {
			return err
		}
		ok4, err := digestsMatch(filepath.Join(calculation, filepath.Base(finalFreq)), freqDigest)
		if err != nil {
			return err
		}
		if !ok1 || !ok2 {
			return fmt.Errorf("Completed dynmat results were modified")
		}
		if !ok3 || !ok4 {
			return fmt.Errorf("Completed dynmat calculation was modified: %s", calculation)
		}
		if err := validateModes(); err != nil {
			return err
		}
		fmt.Println("Preserved completed dynmat.x result")
		return nil
	}
	if _, err := os.Stat(finalOutput); err == nil {
		if isRegularFile(finalOutput) && isRegularFile(finalFreq) {
			for _, calculation := range state.NumberedDirectories(calculations, "dynmat") {
				err := func() error {
					if err := qeoutput.ValidateQeOutput(finalOutput); err != nil {
						return err
					}
					if err := validateModes(); err != nil {
						return err
					}
					if err := validateCalculation(calculation); err != nil {
						return err
					}
					d1, err := state.FileDigest(finalOutput)
					if err != nil {
						return err
					}
					d2, err := state.FileDigest(finalFreq)
					if err != nil {
						return err
					}
					d3, err := state.FileDigest(filepath.Join(calculation, "dynmat.out"))
					if err != nil {
						return err
					}
					d4, err := state.FileDigest(filepath.Join(calculation, filepath.Base(finalFreq)))
					if err != nil {
						return err
					}
					if d1 != d3 || d2 != d4 {
						return fmt.Errorf("digest mismatch")
					}
					if err := config.WriteText(stateMarker, filepath.Base(calculation)+" "+d1+" "+d2+"\n"); err != nil {
						return err
					}
					fmt.Printf("Recovered completed dynmat.x result from %s\n", filepath.Base(calculation))
					return nil
				}()
				if err == nil {
					return nil
				}
				// Try another retained dynmat calculation.
			}
		}
		failed, err := state.NewNumberedDirectory(filepath.Join(s.Workdir, "failed"), "dynmat_publish")
		if err != nil {
			return err
		}
		if _, err := os.Stat(finalOutput); err == nil {
			if err := syscall.Rename(finalOutput, filepath.Join(failed, filepath.Base(finalOutput))); err != nil {
				return fmt.Errorf("Cannot move %s", finalOutput)
			}
		}
		if _, err := os.Stat(finalFreq); err == nil {
			if err := syscall.Rename(finalFreq, filepath.Join(failed, filepath.Base(finalFreq))); err != nil {
				return fmt.Errorf("Cannot move %s", finalFreq)
			}
		}
		fmt.Printf("Preserved incomplete dynmat results in %s\n", config.DisplayPath(failed))
	}
	for _, calculation := range state.NumberedDirectories(calculations, "dynmat") {
		err := func() error {
			if err := validateCalculation(calculation); err != nil {
				return err
			}
			if err := copyFile(filepath.Join(calculation, "dynmat.out"), finalOutput); err != nil {
				return err
			}
			if err := copyFile(filepath.Join(calculation, filepath.Base(finalFreq)), finalFreq); err != nil {
				return err
			}
			d1, err := state.FileDigest(finalOutput)
			if err != nil {
				return err
			}
			d2, err := state.FileDigest(finalFreq)
			if err != nil {
				return err
			}
			if err := config.WriteText(stateMarker, filepath.Base(calculation)+" "+d1+" "+d2+"\n"); err != nil {
				return err
			}
			fmt.Printf("Recovered completed dynmat.x calculation from %s\n", filepath.Base(calculation))
			return nil
		}()
		if err == nil {
			return nil
		}
		// Retain incomplete calculations and create a fresh directory below.
	}
	attempt, err := state.NewNumberedDirectory(calculations, "dynmat")
	if err != nil {
		return err
	}
	if err := copyFile(dyn, filepath.Join(attempt, filepath.Base(dyn))); err != nil {
		return err
	}
	if err := copyFile(filepath.Join(resultsDir, "dynmat.in"), filepath.Join(attempt, "dynmat.in")); err != nil {
		return err
	}
	fmt.Println("Running dynmat.x")
	rc, err := process.ShellRun(s.DynmatCommand+" -inp dynmat.in", attempt, filepath.Join(attempt, "dynmat.out"))
	if err != nil {
		return err
	}
	if rc != 0 {
		return fmt.Errorf("dynmat.x failed with exit code %d", rc)
	}
	if err := qeoutput.ValidateQeOutput(filepath.Join(attempt, "dynmat.out")); err != nil {
		return err
	}
	freq := filepath.Join(attempt, s.OutputPrefix+".freq.out")
	if !isRegularNonempty(freq) {
		return fmt.Errorf("dynmat.x did not produce %s", filepath.Base(freq))
	}
	geometry, err := qevibration.ReadQeDynGeometry(filepath.Join(attempt, filepath.Base(dyn)))
	if err != nil {
		return err
	}
	if _, err := qevibration.ReadQeDynmatModes(freq, len(geometry.Masses)); err != nil {
		return err
	}
	if _, err := os.Stat(finalOutput); err == nil {
		return fmt.Errorf("Refuse to overwrite existing dynmat result")
	}
	if _, err := os.Stat(finalFreq); err == nil {
		return fmt.Errorf("Refuse to overwrite existing dynmat result")
	}
	if err := copyFile(filepath.Join(attempt, "dynmat.out"), finalOutput); err != nil {
		return err
	}
	if err := copyFile(freq, finalFreq); err != nil {
		return err
	}
	d1, err := state.FileDigest(finalOutput)
	if err != nil {
		return err
	}
	d2, err := state.FileDigest(finalFreq)
	if err != nil {
		return err
	}
	return config.WriteText(stateMarker, filepath.Base(attempt)+" "+d1+" "+d2+"\n")
}

func isRegularFile(p string) bool {
	fi, err := os.Stat(p)
	return err == nil && fi.Mode().IsRegular()
}

func isRegularNonempty(p string) bool {
	fi, err := os.Stat(p)
	return err == nil && fi.Mode().IsRegular() && fi.Size() > 0
}

func fileExists(p string) bool {
	_, err := os.Stat(p)
	return err == nil
}

func digestsMatch(p, want string) (bool, error) {
	d, err := state.FileDigest(p)
	if err != nil {
		return false, err
	}
	return d == want, nil
}
