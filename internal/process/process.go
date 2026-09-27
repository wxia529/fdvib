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
	"strings"
	"syscall"
	"time"

	"github.com/wxia529/fdvib/internal/diagnostics"
)

// ShellQuote wraps s in single quotes for POSIX shells, escaping embedded
// quotes, like shell_quote().
func ShellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", "'\\''") + "'"
}

// ShellRun executes cmd with /bin/sh -c in cwd, writing stdout and stderr to
// stdoutPath (truncated), and returns the child exit code.
func ShellRun(cmd, cwd, stdoutPath string) (int, error) {
	return ShellRunWithDiagnostics(cmd, cwd, stdoutPath, nil)
}

// ShellRunWithDiagnostics preserves ShellRun semantics and records shell timing.
func ShellRunWithDiagnostics(cmd, cwd, stdoutPath string, log *diagnostics.Logger) (rc int, resultErr error) {
	end := log.Begin("process", diagnostics.Fields{"command": cmd, "cwd": cwd, "output": stdoutPath})
	defer func() {
		reported := resultErr
		if reported == nil && rc != 0 {
			reported = fmt.Errorf("command exited with code %d", rc)
		}
		end(reported)
	}()
	f, err := os.OpenFile(stdoutPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o644)
	if err != nil {
		log.Event("process.open_failed", diagnostics.Fields{"path": stdoutPath, "error": err.Error()})
		return 126, nil
	}
	defer f.Close()

	sh := exec.Command("/bin/sh", "-c", cmd)
	sh.Dir = cwd
	sh.Stdout = f
	sh.Stderr = f
	err = sh.Start()
	if err == nil {
		started := time.Now()
		log.ProcessStarted(started)
		err = sh.Wait()
		returned := time.Now()
		code := -1
		if sh.ProcessState != nil {
			code = sh.ProcessState.ExitCode()
			if ws, ok := sh.ProcessState.Sys().(syscall.WaitStatus); ok && ws.Signaled() {
				code = 128 + int(ws.Signal())
			}
		}
		log.ProcessReturned(returned, returned.Sub(started), code)
	}
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
