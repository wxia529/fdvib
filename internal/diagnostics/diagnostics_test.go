// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package diagnostics

import (
	"bytes"
	"encoding/json"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestLifecycleAndExclusiveCreation(t *testing.T) {
	path := filepath.Join(t.TempDir(), "debug.jsonl")
	var stderr bytes.Buffer
	log, err := Open(path, "test-version", &stderr)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := Open(path, "test-version", io.Discard); err == nil {
		t.Fatal("existing log was replaced")
	}
	info, err := os.Stat(path)
	if err != nil || info.Mode().Perm() != 0600 {
		t.Fatalf("private permissions: %v %v", info, err)
	}
	log.SetTask("disp_0001_x_p", "disp_0001_x_p_001")
	parent := log.Begin("job", nil)
	if err := log.Do("density.copy", Fields{"bytes": 128}, func() error { return errors.New("copy failed") }); err == nil {
		t.Fatal("lost operation error")
	}
	parent(errors.New("copy failed"))
	log.Close(errors.New("copy failed"))
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var events []Fields
	for _, line := range bytes.Split(bytes.TrimSpace(data), []byte("\n")) {
		var e Fields
		if err := json.Unmarshal(line, &e); err != nil {
			t.Fatal(err)
		}
		events = append(events, e)
	}
	if len(events) != 6 {
		t.Fatalf("events: %s", data)
	}
	for i, e := range events {
		if e["seq"] != float64(i+1) || e["schema_version"] != float64(1) {
			t.Fatalf("ordering: %v", e)
		}
	}
	if events[2]["parent_id"] != events[1]["operation_id"] || events[3]["status"] != "error" || events[3]["job"] != "disp_0001_x_p" {
		t.Fatalf("context: %s", data)
	}
	slowest := events[5]["slowest_operations"].([]any)
	if len(slowest) != 1 || slowest[0].(map[string]any)["phase"] != "density.copy" {
		t.Fatalf("parent double counted: %v", slowest)
	}
	if !strings.Contains(stderr.String(), "density.copy started") {
		t.Fatal("missing start notification")
	}
}

type failingWriter struct{ closed bool }

func (w *failingWriter) Write([]byte) (int, error) { return 0, errors.New("disk full") }
func (w *failingWriter) Close() error              { w.closed = true; return nil }

func TestWriteFailureDoesNotFailOperation(t *testing.T) {
	path := filepath.Join(t.TempDir(), "debug.jsonl")
	var stderr bytes.Buffer
	log, err := Open(path, "test", &stderr)
	if err != nil {
		t.Fatal(err)
	}
	_ = log.file.Close()
	writer := &failingWriter{}
	log.file = writer
	for i := 0; i < 2; i++ {
		if err := log.Do("work", nil, func() error { return nil }); err != nil {
			t.Fatal(err)
		}
	}
	log.Close(nil)
	if !writer.closed || strings.Count(stderr.String(), "log disabled") != 1 {
		t.Fatalf("failure handling: %s", stderr.String())
	}
}

func TestNilLogger(t *testing.T) {
	var log *Logger
	sentinel := errors.New("operation failed")
	if err := log.Do("work", nil, func() error { return sentinel }); err != sentinel {
		t.Fatal("changed error")
	}
	log.Event("ignored", nil)
	log.Close(nil)
}
