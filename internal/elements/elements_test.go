// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package elements

import "testing"

func TestStandardElementSymbol(t *testing.T) {
	cases := []struct {
		label string
		want  string
		ok    bool
	}{
		{"C", "C", true},
		{"C1", "C", true},
		{"O_ads", "O", true},
		{"Fe3", "Fe", true},
		{"H2O", "H", true},
		{"Uu", "U", true},
		{"Xx", "", false},
		{"123", "", false},
	}
	for _, c := range cases {
		got, err := StandardElementSymbol(c.label, "test")
		if c.ok {
			if err != nil || got != c.want {
				t.Errorf("StandardElementSymbol(%q) = %q, %v; want %q", c.label, got, err, c.want)
			}
		} else if err == nil {
			t.Errorf("StandardElementSymbol(%q) should fail", c.label)
		}
	}
}

func TestAtomicNumberFromLabel(t *testing.T) {
	for label, want := range map[string]int{
		"C": 6, "H": 1, "Og": 118, "Na1": 11, "Au": 79,
	} {
		got, err := AtomicNumberFromLabel(label, "test")
		if err != nil || got != want {
			t.Errorf("AtomicNumberFromLabel(%q) = %d, %v; want %d", label, got, err, want)
		}
	}
}
