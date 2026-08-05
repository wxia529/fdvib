// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package run

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"testing"

	"github.com/wxia529/fdvib/internal/settings"
)

const waterScf = `&CONTROL
  calculation = 'scf'
  prefix = 'h2o'
  outdir = '/scratch/qe'
  pseudo_dir = '../pseudo'
  tprnfor = .true.
/
&SYSTEM
  ibrav = 0
  nat = 3
  ntyp = 2
/
&ELECTRONS
  conv_thr = 1.0D-8
/
ATOMIC_SPECIES
  O  15.9994  O.pbe.UPF
  H  1.00794  H.pbe.UPF
ATOMIC_POSITIONS angstrom
  O  0.0  0.0  0.0
  H  0.757  0.586  0.0
  H  -0.757  0.586  0.0
CELL_PARAMETERS angstrom
  10.0  0.0  0.0
  0.0  10.0  0.0
  0.0  0.0  10.0
`

// buildFakeQE compiles the fake pw.x and dynmat.x binaries.
func buildFakeQE(t *testing.T) (pw, dynmat string) {
	t.Helper()
	dir := t.TempDir()
	pw = filepath.Join(dir, "fakeqe-pw")
	dynmat = filepath.Join(dir, "fakeqe-dynmat")
	for _, b := range []struct{ out, pkg string }{
		{pw, "github.com/wxia529/fdvib/internal/fakeqe/pw"},
		{dynmat, "github.com/wxia529/fdvib/internal/fakeqe/dynmat"},
	} {
		cmd := exec.Command("go", "build", "-o", b.out, b.pkg)
		cmd.Env = append(os.Environ(), "CGO_ENABLED=0")
		out, err := cmd.CombinedOutput()
		if err != nil {
			t.Fatalf("go build %s: %v\n%s", b.pkg, err, out)
		}
	}
	return pw, dynmat
}

// makeCase writes scf.in + fdvib.in and returns the settings.
func makeCase(t *testing.T, pw, dynmat string, systemType string) (*settings.Settings, string) {
	t.Helper()
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "scf.in"), []byte(waterScf), 0o644); err != nil {
		t.Fatal(err)
	}
	selected := "1,2,3"
	if systemType == "gas" {
		selected = "all"
	}
	content := fmt.Sprintf("scf_input = scf.in\noutdir = fdvib\nsystem_type = %s\n"+
		"selected_atoms = %s\ndisplacement_angstrom = 0.01\npw_command = %s\n"+
		"prefix = system\nrun_dynmat = true\ndynmat_command = %s\n",
		systemType, selected, pw, dynmat)
	if systemType == "gas" {
		content += "multiplicity = 1\n"
	}
	fdvibIn := filepath.Join(dir, "fdvib.in")
	if err := os.WriteFile(fdvibIn, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	s, err := settings.From(fdvibIn, "")
	if err != nil {
		t.Fatal(err)
	}
	return s, dir
}

// captureStdout runs fn with stdout redirected to a buffer.
func captureStdout(fn func() error) (string, error) {
	old := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		return "", err
	}
	os.Stdout = w
	err = fn()
	w.Close()
	os.Stdout = old
	var buf bytes.Buffer
	_, _ = buf.ReadFrom(r)
	return buf.String(), err
}

// parseDynGMatrix extracts the 3N x 3N matrix block from a dynG file.
func parseDynGMatrix(t *testing.T, path string, nat int) []float64 {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(string(data), "\n")
	start := -1
	for i, line := range lines {
		if strings.Contains(line, "Dynamical  Matrix in cartesian axes") {
			start = i
			break
		}
	}
	if start < 0 {
		t.Fatalf("no matrix block in %s", path)
	}
	n3 := 3 * nat
	h := make([]float64, n3*n3)
	row := 0
	col := 0
	started := false
	for _, line := range lines[start:] {
		fields := strings.Fields(line)
		if len(fields) == 2 && started {
			a, _ := strconv.Atoi(fields[0])
			b, _ := strconv.Atoi(fields[1])
			row = 3 * (a - 1)
			col = 3 * (b - 1)
			continue
		}
		if !started {
			if len(fields) == 2 {
				started = true
				a, _ := strconv.Atoi(fields[0])
				b, _ := strconv.Atoi(fields[1])
				row = 3 * (a - 1)
				col = 3 * (b - 1)
			}
			continue
		}
		if len(fields) >= 6 {
			for k := 0; k < 3; k++ {
				x, err := strconv.ParseFloat(fields[2*k], 64)
				if err != nil {
					t.Fatal(err)
				}
				h[row*n3+col+k] = x
			}
			row++
		}
	}
	return h
}

func TestCalculateFullPipeline(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())

	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatalf("calculate: %v\n%s", err, out)
	}

	workdir := s.Workdir
	// Directory layout.
	for _, p := range []string{
		filepath.Join(workdir, "calculations", "init_scf_001", "scf.in"),
		filepath.Join(workdir, "calculations", "init_scf_001", "scf.out"),
		filepath.Join(workdir, "calculations", "init_scf_001", "out", "h2o.save", "charge-density.dat"),
		filepath.Join(workdir, "calculations", "disp_0001_x_p_001", "pw.in"),
		filepath.Join(workdir, "calculations", "disp_0003_z_m_001", "forces.dat"),
		filepath.Join(workdir, "calculations", "dynmat_001", "dynmat.out"),
		filepath.Join(workdir, "results", "system.dynG"),
		filepath.Join(workdir, "results", "dynmat.in"),
		filepath.Join(workdir, "results", "metadata.dat"),
		filepath.Join(workdir, "results", "dynmat.out"),
		filepath.Join(workdir, "results", "system.freq.out"),
		filepath.Join(workdir, "state", "dataset.state"),
		filepath.Join(workdir, "state", "init_scf.complete"),
		filepath.Join(workdir, "state", "disp_0001_x_p.complete"),
		filepath.Join(workdir, "state", "analyze.complete"),
		filepath.Join(workdir, "state", "dynmat.complete"),
		filepath.Join(workdir, "scf.in.reference"),
		filepath.Join(workdir, "fdvib.in.reference"),
		filepath.Join(workdir, "metadata.dat"),
	} {
		if _, err := os.Stat(p); err != nil {
			t.Errorf("missing %s", p)
		}
	}
	// 18 displacement attempts, one job each.
	entries, _ := os.ReadDir(filepath.Join(workdir, "calculations"))
	dispCount := 0
	for _, e := range entries {
		if strings.HasPrefix(e.Name(), "disp_") {
			dispCount++
		}
	}
	if dispCount != 18 {
		t.Errorf("displacement attempts = %d, want 18", dispCount)
	}
	// Displacement densities must have been removed.
	for _, e := range entries {
		if strings.HasPrefix(e.Name(), "disp_") {
			if _, err := os.Stat(filepath.Join(workdir, "calculations", e.Name(), "out", "h2o.save", "charge-density.dat")); err == nil {
				t.Errorf("displacement density not removed in %s", e.Name())
			}
		}
	}
	// dynG matrix must reproduce the fake Hessian.
	got := parseDynGMatrix(t, filepath.Join(workdir, "results", "system.dynG"), 3)
	want := fakeHessian(3, 0.4, 0.05, 0.02)
	for i := range got {
		if math.Abs(got[i]-want[i]) > 1e-6 {
			t.Errorf("dynG[%d] = %v, want %v", i, got[i], want[i])
		}
	}
	// dynmat.in content.
	di, _ := os.ReadFile(filepath.Join(workdir, "results", "dynmat.in"))
	diText := string(di)
	for _, want := range []string{"fildyn='system.dynG',", "filout='system.freq.out',",
		"asr='no',", "remove_interaction_blocks=.true.,"} {
		if !strings.Contains(diText, want) {
			t.Errorf("dynmat.in missing %q:\n%s", want, diText)
		}
	}
	// metadata.dat has 16-decimal scientific energy.
	md, _ := os.ReadFile(filepath.Join(workdir, "metadata.dat"))
	mdText := string(md)
	if !strings.Contains(mdText, "electronic_energy_hartree = -8.5617283945061704e+00") {
		t.Errorf("metadata energy wrong:\n%s", mdText)
	}
	if !strings.Contains(mdText, "mode_selection = local\nselected_atoms = 1,2,3") {
		t.Errorf("metadata selection wrong:\n%s", mdText)
	}
	// freq.out has 9 modes.
	fo, _ := os.ReadFile(filepath.Join(workdir, "results", "system.freq.out"))
	if strings.Count(string(fo), "freq (") != 9 {
		t.Errorf("freq.out modes:\n%s", fo)
	}
	// dataset.state format.
	ds, _ := os.ReadFile(filepath.Join(workdir, "state", "dataset.state"))
	if !strings.HasPrefix(string(ds), "format=2\n") {
		t.Errorf("dataset.state:\n%s", ds)
	}
}

func fakeHessian(nat int, a, b, c float64) []float64 {
	n3 := 3 * nat
	h := make([]float64, n3*n3)
	for i := 0; i < nat; i++ {
		for j := 0; j < nat; j++ {
			for alpha := 0; alpha < 3; alpha++ {
				for beta := 0; beta < 3; beta++ {
					row := 3*i + alpha
					col := 3*j + beta
					var v float64
					if i == j {
						v = b
						if alpha == beta {
							v += a
						}
					} else if abs(i-j) == 1 {
						v = c
					}
					h[row*n3+col] = v
				}
			}
		}
	}
	return h
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func TestCalculateRerunPreserves(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	if _, err := captureStdout(func() error { return Calculate(s) }); err != nil {
		t.Fatal(err)
	}
	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{
		"Preserved completed reference SCF",
		"preserved 18 displacement jobs",
		"Preserved completed Hessian analysis",
		"Preserved completed dynmat.x result",
		"FDVIB calculation completed",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("rerun output missing %q:\n%s", want, out)
		}
	}
	// No new attempt directories.
	entries, _ := os.ReadDir(filepath.Join(s.Workdir, "calculations"))
	if len(entries) != 20 { // init_scf + 18 disp + dynmat
		t.Errorf("attempts = %d", len(entries))
	}
}

func TestCalculateRecovery(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	if _, err := captureStdout(func() error { return Calculate(s) }); err != nil {
		t.Fatal(err)
	}
	workdir := s.Workdir
	// Remove the analyze marker: results exist, so it must be recovered.
	if err := os.Remove(filepath.Join(workdir, "state", "analyze.complete")); err != nil {
		t.Fatal(err)
	}
	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out, "Recovered completed Hessian analysis") {
		t.Errorf("missing recovery message:\n%s", out)
	}
	if _, err := os.Stat(filepath.Join(workdir, "state", "analyze.complete")); err != nil {
		t.Errorf("analyze marker not rewritten: %v", err)
	}

	// Remove one displacement marker: the job must be recovered from its
	// retained directory.
	if err := os.Remove(filepath.Join(workdir, "state", "disp_0001_x_p.complete")); err != nil {
		t.Fatal(err)
	}
	out, err = captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out, "Recovered completed disp_0001_x_p from disp_0001_x_p_001") {
		t.Errorf("missing displacement recovery:\n%s", out)
	}
}

func TestCalculateQuarantinesBrokenResults(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	if _, err := captureStdout(func() error { return Calculate(s) }); err != nil {
		t.Fatal(err)
	}
	workdir := s.Workdir
	// Corrupt the dynG and drop the analyze marker: analysis must rebuild
	// from scratch, moving the old results to failed/. The dynmat marker is
	// dropped too because the rebuilt results no longer contain dynmat.out.
	if err := os.Remove(filepath.Join(workdir, "state", "analyze.complete")); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove(filepath.Join(workdir, "state", "dynmat.complete")); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(workdir, "results", "system.dynG"), []byte("garbage\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatalf("calculate: %v\n%s", err, out)
	}
	if !strings.Contains(out, "Preserved incomplete Hessian results in") {
		t.Errorf("missing quarantine message:\n%s", out)
	}
	entries, _ := os.ReadDir(filepath.Join(workdir, "failed"))
	found := false
	for _, e := range entries {
		if strings.HasPrefix(e.Name(), "analysis_") {
			found = true
		}
	}
	if !found {
		t.Error("failed/analysis_* missing")
	}
	// The rebuilt dynG must again match the fake Hessian.
	got := parseDynGMatrix(t, filepath.Join(workdir, "results", "system.dynG"), 3)
	want := fakeHessian(3, 0.4, 0.05, 0.02)
	for i := range got {
		if math.Abs(got[i]-want[i]) > 1e-6 {
			t.Errorf("dynG[%d] = %v, want %v", i, got[i], want[i])
		}
	}
}

func TestCalculateQuarantinesBrokenDynmat(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	if _, err := captureStdout(func() error { return Calculate(s) }); err != nil {
		t.Fatal(err)
	}
	workdir := s.Workdir
	// Corrupt the published freq.out and drop the dynmat marker: the stale
	// publish is quarantined and dynmat.x reruns.
	if err := os.Remove(filepath.Join(workdir, "state", "dynmat.complete")); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(workdir, "results", "system.freq.out"), []byte("junk\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatalf("calculate: %v\n%s", err, out)
	}
	if !strings.Contains(out, "Preserved incomplete dynmat results in") {
		t.Errorf("missing quarantine message:\n%s", out)
	}
	if !strings.Contains(out, "Recovered completed dynmat.x calculation from dynmat_001") {
		t.Errorf("dynmat.x calculation not recovered:\n%s", out)
	}
	fo, err := os.ReadFile(filepath.Join(workdir, "results", "system.freq.out"))
	if err != nil || !strings.Contains(string(fo), "freq (") {
		t.Errorf("freq.out not republished: %v", err)
	}
}

func TestCalculateLock(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	// Take the lock ourselves; Calculate must refuse.
	lockPath := s.Workdir + ".lock"
	f, err := os.OpenFile(lockPath, os.O_CREATE|os.O_RDWR, 0o644)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	if err := syscall.Flock(int(f.Fd()), syscall.LOCK_EX|syscall.LOCK_NB); err != nil {
		t.Fatal(err)
	}
	err = Calculate(s)
	if err == nil || !strings.Contains(err.Error(), "Calculation is already running") {
		t.Errorf("expected lock refusal, got %v", err)
	}
	syscall.Flock(int(f.Fd()), syscall.LOCK_UN)
	if _, err := captureStdout(func() error { return Calculate(s) }); err != nil {
		t.Fatalf("after unlock: %v", err)
	}
}

func TestCalculateRunDynmatFalse(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "scf.in"), []byte(waterScf), 0o644); err != nil {
		t.Fatal(err)
	}
	content := fmt.Sprintf("scf_input = scf.in\noutdir = fdvib\nsystem_type = local\n"+
		"selected_atoms = 1,2,3\ndisplacement_angstrom = 0.01\npw_command = %s\n"+
		"prefix = system\nrun_dynmat = false\ndynmat_command = %s\n", pw, dynmat)
	fdvibIn := filepath.Join(dir, "fdvib.in")
	if err := os.WriteFile(fdvibIn, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	s, err := settings.From(fdvibIn, "")
	if err != nil {
		t.Fatal(err)
	}
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatalf("calculate: %v", err)
	}
	if !strings.Contains(out, "dynmat.x was not requested (run_dynmat=.false.)") {
		t.Errorf("missing message:\n%s", out)
	}
	if _, err := os.Stat(filepath.Join(s.Workdir, "results", "system.freq.out")); err == nil {
		t.Error("freq.out should not exist")
	}
	// Rerun with run_dynmat=true: dynmat.x runs without redoing SCF work.
	content = strings.Replace(content, "run_dynmat = false", "run_dynmat = true", 1)
	if err := os.WriteFile(fdvibIn, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	s2, err := settings.From(fdvibIn, "")
	if err != nil {
		t.Fatal(err)
	}
	out, err = captureStdout(func() error { return Calculate(s2) })
	if err != nil {
		t.Fatalf("calculate: %v", err)
	}
	if !strings.Contains(out, "Running dynmat.x") {
		t.Errorf("dynmat.x not run:\n%s", out)
	}
	if _, err := os.Stat(filepath.Join(s2.Workdir, "results", "system.freq.out")); err != nil {
		t.Error("freq.out missing after enabling dynmat")
	}
}

func TestCalculateGas(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "gas")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	out, err := captureStdout(func() error { return Calculate(s) })
	if err != nil {
		t.Fatalf("calculate: %v\n%s", err, out)
	}
	if !strings.Contains(out, "FDVIB calculation completed") {
		t.Errorf("output:\n%s", out)
	}
	md, _ := os.ReadFile(filepath.Join(s.Workdir, "results", "metadata.dat"))
	if !strings.Contains(string(md), "mode_selection = gas\nselected_atoms = all") {
		t.Errorf("metadata:\n%s", md)
	}
	di, _ := os.ReadFile(filepath.Join(s.Workdir, "results", "dynmat.in"))
	if !strings.Contains(string(di), "remove_interaction_blocks=.false.,") {
		t.Errorf("dynmat.in:\n%s", di)
	}
}

func TestCalculateGasValidation(t *testing.T) {
	pw, _ := buildFakeQE(t)
	dir := t.TempDir()
	os.WriteFile(filepath.Join(dir, "scf.in"), []byte(waterScf), 0o644)
	content := fmt.Sprintf("scf_input = scf.in\noutdir = fdvib\nsystem_type = gas\n"+
		"selected_atoms = all\ndisplacement_angstrom = 0.01\npw_command = %s\n"+
		"prefix = system\nrun_dynmat = false\n", pw)
	// No multiplicity: gas requires it explicitly.
	fdvibIn := filepath.Join(dir, "fdvib.in")
	os.WriteFile(fdvibIn, []byte(content), 0o644)
	s, err := settings.From(fdvibIn, "")
	if err != nil {
		t.Fatal(err)
	}
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	if err := Calculate(s); err == nil || !strings.Contains(err.Error(), "gas requires an explicit positive multiplicity") {
		t.Errorf("expected multiplicity error, got %v", err)
	}
	// Multiplicity=2 but the QE input has no nspin/tot_magnetization.
	content += "multiplicity = 2\n"
	os.WriteFile(fdvibIn, []byte(content), 0o644)
	s, err = settings.From(fdvibIn, "")
	if err != nil {
		t.Fatal(err)
	}
	if err := Calculate(s); err == nil || !strings.Contains(err.Error(), "nspin=2 and tot_magnetization=multiplicity-1") {
		t.Errorf("expected spin error, got %v", err)
	}
}

func TestCalculateRefusesChangedDataset(t *testing.T) {
	pw, dynmat := buildFakeQE(t)
	s, _ := makeCase(t, pw, dynmat, "local")
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	if _, err := captureStdout(func() error { return Calculate(s) }); err != nil {
		t.Fatal(err)
	}
	// Modify scf.in: dataset differs.
	scf := filepath.Join(s.Root, "scf.in")
	data, _ := os.ReadFile(scf)
	if err := os.WriteFile(scf, append(data, []byte("\n! changed\n")...), 0o644); err != nil {
		t.Fatal(err)
	}
	err := Calculate(s)
	if err == nil || !strings.Contains(err.Error(), "Dataset differs from the existing calculation") {
		t.Errorf("expected dataset error, got %v", err)
	}
}

func TestCalculateRefusesNonEmptyOutdirWithoutState(t *testing.T) {
	pw, _ := buildFakeQE(t)
	s, _ := makeCase(t, pw, "/nonexistent/dynmat.x", "local")
	if err := os.MkdirAll(filepath.Join(s.Workdir, "junk"), 0o755); err != nil {
		t.Fatal(err)
	}
	os.WriteFile(filepath.Join(s.Workdir, "junk", "x"), []byte("x"), 0o644)
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	err := Calculate(s)
	if err == nil || !strings.Contains(err.Error(), "Refusing to use non-empty outdir without FDVIB state metadata") {
		t.Errorf("expected refusal, got %v", err)
	}
}

func TestCalculateFailureExitCode(t *testing.T) {
	// pw_command that fails: the reference SCF reports the exit code.
	dir := t.TempDir()
	os.WriteFile(filepath.Join(dir, "scf.in"), []byte(waterScf), 0o644)
	content := "scf_input = scf.in\noutdir = fdvib\nsystem_type = local\n" +
		"selected_atoms = 1,2,3\ndisplacement_angstrom = 0.01\n" +
		"pw_command = /bin/false\nprefix = system\nrun_dynmat = false\n"
	fdvibIn := filepath.Join(dir, "fdvib.in")
	os.WriteFile(fdvibIn, []byte(content), 0o644)
	s, err := settings.From(fdvibIn, "")
	if err != nil {
		t.Fatal(err)
	}
	t.Setenv("FDVIB_FAKE_STATE", t.TempDir())
	err = Calculate(s)
	if err == nil || !strings.Contains(err.Error(), "Reference SCF failed with exit code 1") {
		t.Errorf("expected exit-code error, got %v", err)
	}
}
