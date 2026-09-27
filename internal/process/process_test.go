// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package process

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/wxia529/fdvib/internal/diagnostics"
)

func TestWaitsAfterJobDone(t *testing.T) {
	dir := t.TempDir()
	output := filepath.Join(dir, "pw.out")
	release := filepath.Join(dir, "release")
	logPath := filepath.Join(dir, "debug.jsonl")
	log, err := diagnostics.Open(logPath, "test", io.Discard)
	if err != nil {
		t.Fatal(err)
	}
	defer log.Close(nil)
	done := make(chan error, 1)
	finished := make(chan struct{})
	go func() {
		defer close(finished)
		rc, err := ShellRunWithDiagnostics("printf 'JOB DONE\\n'; while [ ! -f release ]; do sleep 0.01; done", dir, output, log)
		if err == nil && rc != 0 {
			err = io.ErrUnexpectedEOF
		}
		done <- err
	}()
	// Always release the shell, including when the assertion fails.
	defer func() {
		_ = os.WriteFile(release, nil, 0600)
		select {
		case <-finished:
		case <-time.After(5 * time.Second):
		}
	}()
	deadline := time.After(5 * time.Second)
	for {
		data, _ := os.ReadFile(output)
		if bytes.Contains(data, []byte("JOB DONE")) {
			break
		}
		select {
		case <-deadline:
			t.Fatal("shell did not produce JOB DONE")
		case <-time.After(10 * time.Millisecond):
		}
	}
	select {
	case err := <-done:
		t.Fatalf("returned before shell exited: %v", err)
	default:
	}
	if err := os.WriteFile(release, nil, 0600); err != nil {
		t.Fatal(err)
	}
	select {
	case err := <-done:
		if err != nil {
			t.Fatal(err)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("shell did not exit")
	}
	log.Close(nil)
	data, err := os.ReadFile(logPath)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Contains(data, []byte(`"event":"process.return"`)) {
		t.Fatalf("missing return: %s", data)
	}
}

func TestExitSemantics(t *testing.T) {
	for _, tc := range []struct {
		command string
		code    int
	}{{"exit 7", 7}, {"exec /definitely/missing/program", 127}, {"kill -TERM $$", 143}} {
		dir := t.TempDir()
		rc, err := ShellRun(tc.command, dir, filepath.Join(dir, "out"))
		if err != nil || rc != tc.code {
			t.Fatalf("%s: code=%d error=%v", tc.command, rc, err)
		}
	}
	dir := t.TempDir()
	rc, err := ShellRun("true", dir, filepath.Join(dir, "missing", "out"))
	if rc != 126 || err != nil {
		t.Fatalf("open failure: %d %v", rc, err)
	}
}
