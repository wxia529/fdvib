// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package diagnostics provides optional, local operation timing logs.
package diagnostics

import (
	"encoding/json"
	"fmt"
	"io"
	"os"
	"runtime"
	"sort"
	"time"
)

type Fields map[string]any

type sample struct {
	Phase      string  `json:"phase"`
	DurationMS float64 `json:"duration_ms"`
	Job        string  `json:"job,omitempty"`
	Attempt    string  `json:"attempt,omitempty"`
}

type total struct {
	Count      int     `json:"count"`
	DurationMS float64 `json:"duration_ms"`
}

// Logger is scoped to one invocation. A nil Logger disables diagnostics.
// Operations are synchronous, matching the calculation orchestrator.
type Logger struct {
	job, attempt         string
	file                 io.WriteCloser
	stderr               io.Writer
	start                time.Time
	runID                string
	seq, next            uint64
	stack                []uint64
	totals               map[string]total
	slowest              []sample
	external             time.Duration
	lastReturn           time.Time
	lastJob, lastAttempt string
	closed               bool
}

// Open creates a private log without replacing existing files or creating
// parent directories. Empty path selects a unique file in the current directory.
func Open(path, version string, stderr io.Writer) (*Logger, error) {
	var f *os.File
	var err error
	if path == "" {
		f, err = os.CreateTemp(".", "fdvib-debug-*.jsonl")
	} else {
		f, err = os.OpenFile(path, os.O_WRONLY|os.O_CREATE|os.O_EXCL, 0600)
	}
	if err != nil {
		return nil, fmt.Errorf("cannot create debug log: %w", err)
	}
	now := time.Now()
	l := &Logger{file: f, stderr: stderr, start: now, runID: fmt.Sprintf("%d-%d", os.Getpid(), now.UnixNano()), totals: make(map[string]total)}
	fmt.Fprintf(stderr, "[fdvib debug] log: %s\n", f.Name())
	cwd, _ := os.Getwd()
	l.Event("run.start", Fields{"version": version, "go_version": runtime.Version(), "pid": os.Getpid(), "cwd": cwd})
	return l, nil
}

func (l *Logger) write(fields Fields) {
	if l.file == nil {
		return
	}
	l.seq++
	if l.job != "" {
		if _, ok := fields["job"]; !ok {
			fields["job"] = l.job
		}
	}
	if l.attempt != "" {
		if _, ok := fields["attempt"]; !ok {
			fields["attempt"] = l.attempt
		}
	}
	fields["schema_version"], fields["run_id"], fields["seq"] = 1, l.runID, l.seq
	fields["time"] = time.Now().UTC().Format(time.RFC3339Nano)
	fields["elapsed_ms"] = float64(time.Since(l.start)) / float64(time.Millisecond)
	data, err := json.Marshal(fields)
	if err == nil {
		data = append(data, '\n')
		var n int
		n, err = l.file.Write(data)
		if err == nil && n != len(data) {
			err = io.ErrShortWrite
		}
	}
	if err != nil {
		fmt.Fprintf(l.stderr, "[fdvib debug] log disabled after write failure: %v\n", err)
		_ = l.file.Close()
		l.file = nil
	}
}

func clone(fields Fields) Fields {
	result := make(Fields, len(fields)+8)
	for k, v := range fields {
		result[k] = v
	}
	return result
}

// SetTask adds task context to subsequent events. Empty values clear it.
func (l *Logger) SetTask(job, attempt string) {
	if l != nil {
		l.job, l.attempt = job, attempt
	}
}

// Event records a point event without timing it.
func (l *Logger) Event(event string, fields Fields) {
	if l == nil {
		return
	}
	f := clone(fields)
	f["event"] = event
	if len(l.stack) > 0 {
		f["parent_id"] = l.stack[len(l.stack)-1]
	}
	l.write(f)
}

// Begin records an operation immediately. The returned function must be
// called exactly once, in nesting order, with the operation's final error.
func (l *Logger) Begin(phase string, fields Fields) func(error) {
	if l == nil {
		return func(error) {}
	}
	l.next++
	id := l.next
	f := clone(fields)
	f["event"], f["phase"], f["operation_id"] = "operation.start", phase, id
	if len(l.stack) > 0 {
		f["parent_id"] = l.stack[len(l.stack)-1]
	}
	start := time.Now()
	l.write(f)
	label := phase
	if l.job != "" {
		label = l.job + ": " + phase
	}
	fmt.Fprintf(l.stderr, "[fdvib debug] %s started\n", label)
	l.stack = append(l.stack, id)
	return func(err error) {
		duration := time.Since(start)
		l.stack = l.stack[:len(l.stack)-1]
		end := clone(f)
		end["event"], end["duration_ms"], end["status"] = "operation.end", float64(duration)/float64(time.Millisecond), "success"
		if err != nil {
			end["status"], end["error"] = "error", err.Error()
		}
		l.write(end)
		fmt.Fprintf(l.stderr, "[fdvib debug] %s %s, %.3f s\n", label, end["status"], duration.Seconds())
		t := l.totals[phase]
		t.Count++
		t.DurationMS += float64(duration) / float64(time.Millisecond)
		l.totals[phase] = t
		if id == l.next {
			l.slowest = append(l.slowest, sample{Phase: phase, DurationMS: float64(duration) / float64(time.Millisecond), Job: l.job, Attempt: l.attempt})
		}
		sort.Slice(l.slowest, func(i, j int) bool { return l.slowest[i].DurationMS > l.slowest[j].DurationMS })
		if len(l.slowest) > 5 {
			l.slowest = l.slowest[:5]
		}
	}
}

// Do times a synchronous operation.
func (l *Logger) Do(phase string, fields Fields, fn func() error) error {
	end := l.Begin(phase, fields)
	err := fn()
	end(err)
	return err
}

// ProcessStarted and ProcessReturned measure shell lifetime and the gap
// between successive shell processes, not QE's JOB DONE output marker.
func (l *Logger) ProcessStarted(started time.Time) {
	if l == nil {
		return
	}
	l.Event("process.started", Fields{"started_time": started.UTC().Format(time.RFC3339Nano)})
	if !l.lastReturn.IsZero() {
		l.Event("process.transition", Fields{"duration_ms": float64(started.Sub(l.lastReturn)) / float64(time.Millisecond), "from_job": l.lastJob, "from_attempt": l.lastAttempt})
	}
}

func (l *Logger) ProcessReturned(returned time.Time, duration time.Duration, rc int) {
	if l == nil {
		return
	}
	l.external += duration
	l.lastReturn = returned
	l.lastJob, l.lastAttempt = l.job, l.attempt
	l.Event("process.return", Fields{"duration_ms": float64(duration) / float64(time.Millisecond), "exit_code": rc})
}

// Close records the final status. Logging failures never change calculation success.
func (l *Logger) Close(err error) {
	if l == nil || l.closed {
		return
	}
	l.closed = true
	status := "success"
	if err != nil {
		status = "error"
	}
	elapsed := time.Since(l.start)
	fields := Fields{"status": status, "duration_ms": float64(elapsed) / float64(time.Millisecond), "external_process_ms": float64(l.external) / float64(time.Millisecond), "fdvib_ms": float64(elapsed-l.external) / float64(time.Millisecond), "phases": l.totals, "slowest_operations": l.slowest}
	if err != nil {
		fields["error"] = err.Error()
	}
	l.Event("run.end", fields)
	fmt.Fprintf(l.stderr, "[fdvib debug] run %s, total %.3f s, external commands %.3f s, FDVIB %.3f s\n", status, elapsed.Seconds(), l.external.Seconds(), (elapsed - l.external).Seconds())
	if l.file != nil {
		if closeErr := l.file.Close(); closeErr != nil {
			fmt.Fprintf(l.stderr, "[fdvib debug] cannot close log: %v\n", closeErr)
		}
		l.file = nil
	}
}
