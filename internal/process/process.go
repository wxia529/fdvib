// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package process runs QE commands through /bin/sh with the same exit-code
// semantics as shell_run() in process.cpp (126 chdir/open failure, 127 exec
// failure, 128+signal), capturing stdout and stderr to a file.
package process

import (
	"errors"
	"fmt"
	"os"
	"os/exec"
	"syscall"
)

// ShellRun executes cmd with /bin/sh -c in cwd, writing stdout and stderr to
// stdoutPath (truncated), and returns the child exit code.
func ShellRun(cmd, cwd, stdoutPath string) (int, error) {
	f, err := os.OpenFile(stdoutPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o644)
	if err != nil {
		return 126, nil
	}
	defer f.Close()

	sh := exec.Command("/bin/sh", "-c", cmd)
	sh.Dir = cwd
	sh.Stdout = f
	sh.Stderr = f
	err = sh.Run()
	if err == nil {
		return 0, nil
	}
	if ee, ok := err.(*exec.ExitError); ok {
		if ws, ok := ee.Sys().(syscall.WaitStatus); ok {
			if ws.Signaled() {
				return 128 + int(ws.Signal()), nil
			}
			return ws.ExitStatus(), nil
		}
		return ee.ExitCode(), nil
	}
	if errors.Is(err, exec.ErrNotFound) {
		return 127, nil
	}
	// Process could not start (e.g. cwd missing), matching the child's 126.
	return 126, fmt.Errorf("cannot start command: %w", err)
}
