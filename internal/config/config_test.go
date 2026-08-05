// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

package config

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
	"testing"
)

// TestFormatMatchesCPP verifies that FormatGeneral/FormatSci/FormatFixed
// match the frozen formatting baseline byte-for-byte for a broad set of
// values and precisions; fmt_cpp.txt is the frozen format baseline.
func TestFormatMatchesCPP(t *testing.T) {
	f, err := os.Open("testdata/fmt_cpp.txt")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	values := []float64{
		0.0, math.Copysign(0, -1), 1.0, -1.0, 1.5, 0.5, 0.001, 123.45, 999.9999999999999,
		1e5, 1e-5, 1.23e-7, -3.5, 123456789.0, 0.1, 0.01, 1.0 / 3.0, 2.5e10,
		3e-10, 1e15, 123456789012345.0, 0.000123456789, 6.02214076e23,
		5.0e-324, 2.2250738585072014e-308, 1.7976931348623157e308, 1000.0,
		1000000000000000.0, 123456789.123456789, 0.0001, 0.00009999999999999999,
	}
	precs := []int{6, 7, 9, 10, 12, 17}

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)
	lineNo := 0
	for _, x := range values {
		for _, p := range precs {
			want := scanLine(t, scanner, &lineNo)
			got := "G" + strconv.Itoa(p) + ":" + FormatGeneral(x, p)
			if got != want {
				t.Errorf("FormatGeneral(%v, %d) = %q, want %q", x, p, got, want)
			}
		}
		for _, p := range precs {
			want := scanLine(t, scanner, &lineNo)
			got := "S" + strconv.Itoa(p) + ":" + FormatSci(x, p)
			if got != want {
				t.Errorf("FormatSci(%v, %d) = %q, want %q", x, p, got, want)
			}
		}
		for _, p := range precs {
			want := scanLine(t, scanner, &lineNo)
			got := "F" + strconv.Itoa(p) + ":" + FormatFixed(x, p)
			if got != want {
				t.Errorf("FormatFixed(%v, %d) = %q, want %q", x, p, got, want)
			}
		}
	}
	if scanner.Scan() {
		t.Fatalf("reference file has extra lines: %q", scanner.Text())
	}
}

func scanLine(t *testing.T, scanner *bufio.Scanner, lineNo *int) string {
	t.Helper()
	if !scanner.Scan() {
		t.Fatalf("reference file ended early at line %d", *lineNo)
	}
	*lineNo++
	return strings.TrimSuffix(scanner.Text(), "\n")
}

func TestTrimLowerUnquote(t *testing.T) {
	if got := strings.Trim("  a b \t\r\n", " \t\r\n"); got != "a b" {
		t.Errorf("Trim = %q", got)
	}
	if got := strings.ToLower("AbC"); got != "abc" {
		t.Errorf("Lower = %q", got)
	}
	for in, want := range map[string]string{
		"'x'": "x", `"y"`: "y", " z ": "z", "'a'b'": "a'b",
		"nope": "nope", "''": "", `""`: "",
	} {
		if got := Unquote(in); got != want {
			t.Errorf("Unquote(%q) = %q, want %q", in, got, want)
		}
	}
}

func TestStripComment(t *testing.T) {
	// C++ strip_comment cuts at '!' only (outside quotes); '#' is only a
	// comment marker in the key-value parser (strip_key_value_comment).
	cases := []struct{ in, want string }{
		{"a ! b", "a "},
		{"'!' ! x", "'!' "},
		{`"!" ! x`, `"!" `},
		{"no comment", "no comment"},
		{"a # b", "a # b"},
	}
	for _, c := range cases {
		if got := StripComment(c.in); got != c.want {
			t.Errorf("StripComment(%q) = %q, want %q", c.in, got, c.want)
		}
	}
}

func TestParseNumber(t *testing.T) {
	for in, want := range map[string]float64{
		"1.5": 1.5, "1D-3": 1e-3, "1d+2": 100.0, "-2.5e1": -25.0,
		"'3.0'": 3.0, " .5 ": 0.5, "1e2": 100.0,
	} {
		got, err := ParseNumber(in)
		if err != nil || got != want {
			t.Errorf("ParseNumber(%q) = %v, %v; want %v", in, got, err, want)
		}
	}
	for _, in := range []string{"1.2.3", "abc", "inf", "nan", "1e999"} {
		if _, err := ParseNumber(in); err == nil {
			t.Errorf("ParseNumber(%q) should fail", in)
		}
	}
}

func TestParseIntList(t *testing.T) {
	got, err := ParseIntList("1,2;3 4")
	if err != nil || len(got) != 4 || got[0] != 1 || got[3] != 4 {
		t.Errorf("IntegerList = %v, %v", got, err)
	}
	if got, err := ParseIntList(""); err != nil || len(got) != 0 {
		t.Errorf("ParseIntList('') = %v, %v", got, err)
	}
	if _, err := ParseIntList("1 x"); err == nil {
		t.Error("IntegerList should reject non-integers")
	}
}

func TestFixedField(t *testing.T) {
	// FixedField = right-padded fixed formatting with width, like
	// `setw(width) << fixed << setprecision(prec) << v`.
	cases := []struct {
		v           float64
		prec, width int
		want        string
	}{
		{1.5, 2, 8, "    1.50"},
		{1.5, 2, 4, "1.50"},
		{0, 7, 12, "   0.0000000"},
		{-3.25, 3, 10, "    -3.250"},
		{math.Inf(1), 2, 6, "   inf"},
		{math.NaN(), 2, 6, "   nan"},
	}
	for _, c := range cases {
		if got := FixedField(c.v, c.prec, c.width); got != c.want {
			t.Errorf("FixedField(%v, %d, %d) = %q, want %q", c.v, c.prec, c.width, got, c.want)
		}
	}
}
