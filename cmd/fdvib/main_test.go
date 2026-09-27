// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package main

import (
	"reflect"
	"testing"
)

func TestDebugArgs(t *testing.T) {
	for _, args := range [][]string{
		{"fdvib", "--debug", "-inp", "fdvib.in", "--debug-log", "debug.jsonl"},
		{"fdvib", "-inp", "fdvib.in", "--debug-log", "debug.jsonl", "--debug"},
	} {
		clean, enabled, path, err := debugArgs(args)
		if err != nil || !enabled || path != "debug.jsonl" || !reflect.DeepEqual(clean, []string{"fdvib", "-inp", "fdvib.in"}) {
			t.Fatalf("%v: %v %t %s %v", args, clean, enabled, path, err)
		}
	}
	for _, args := range [][]string{
		{"fdvib", "--debug-log", "x"},
		{"fdvib", "--debug", "--debug"},
		{"fdvib", "--debug", "--debug-log"},
		{"fdvib", "--debug", "--debug-log", "--help"},
		{"fdvib", "--debug", "--debug-log", "a", "--debug-log", "b"},
	} {
		if _, _, _, err := debugArgs(args); err == nil {
			t.Fatalf("accepted %v", args)
		}
	}
}
