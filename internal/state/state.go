// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package state provides FNV-1a file digests and the numbered attempt
// directory scheme used for recovery, mirroring run_state.cpp.
package state

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// fnv1a64 implements the standard FNV-1a 64-bit hash, matching the C++
// implementation byte-for-byte. Go's hash/fnv package uses a different
// (non-standard) offset basis and cannot be used here.
func fnv1a64(data []byte) uint64 {
	const (
		offset64 = uint64(1469598103934665603)
		prime64  = uint64(1099511628211)
	)
	var hash = offset64
	for _, b := range data {
		hash ^= uint64(b)
		hash *= prime64
	}
	return hash
}

// FileDigest computes the FNV-1a 64-bit digest of a file as 16 lowercase hex
// digits, like file_digest().
func FileDigest(p string) (string, error) {
	data, err := os.ReadFile(p)
	if err != nil {
		return "", fmt.Errorf("cannot read %s", p)
	}
	return fmt.Sprintf("%016x", fnv1a64(data)), nil
}

// NumberedName formats a task attempt name like task_001, like
// numbered_name().
func NumberedName(task string, n int) string {
	return fmt.Sprintf("%s_%03d", task, n)
}

// IsNumberedName reports whether name is task_NNN with at least three
// digits, like is_numbered_name().
func IsNumberedName(name, task string) bool {
	prefix := task + "_"
	if !strings.HasPrefix(name, prefix) || len(name) < len(prefix)+3 {
		return false
	}
	for _, c := range name[len(prefix):] {
		if c < '0' || c > '9' {
			return false
		}
	}
	return true
}

// NewNumberedDirectory creates parent and the next free task_NNN directory,
// like new_numbered_directory().
func NewNumberedDirectory(parent, task string) (string, error) {
	if err := os.MkdirAll(parent, 0o755); err != nil {
		return "", fmt.Errorf("cannot create %s", parent)
	}
	n := 1
	for {
		name := NumberedName(task, n)
		path := filepath.Join(parent, name)
		if _, err := os.Stat(path); err != nil {
			if err := os.MkdirAll(path, 0o755); err != nil {
				return "", fmt.Errorf("cannot create %s", path)
			}
			return path, nil
		}
		n++
	}
}

// NumberedDirectories lists the task_NNN directories under parent sorted by
// name descending, like numbered_directories().
func NumberedDirectories(parent, task string) []string {
	var found []string
	info, err := os.Stat(parent)
	if err != nil || !info.IsDir() {
		return found
	}
	entries, err := os.ReadDir(parent)
	if err != nil {
		return found
	}
	for _, entry := range entries {
		if entry.IsDir() && IsNumberedName(entry.Name(), task) {
			found = append(found, filepath.Join(parent, entry.Name()))
		}
	}
	// Descending string order, like std::sort(rbegin, rend).
	for i := 0; i < len(found); i++ {
		for j := i + 1; j < len(found); j++ {
			if found[j] > found[i] {
				found[i], found[j] = found[j], found[i]
			}
		}
	}
	return found
}
