// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package config provides text utilities, the strict key=value configuration
// parser, and the numeric formatting helpers that reproduce C++ iostream
// output semantics byte-for-byte. It mirrors common.cpp and the constants in
// fdvib.hpp.
package config

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// Trim strips whitespace (" \t\r\n") from both ends of s, like C++ trim().
func Trim(s string) string {
	return strings.Trim(s, " \t\r\n")
}

// Lower lowercases s, like C++ lower().
func Lower(s string) string {
	return strings.ToLower(s)
}

// Unquote removes one layer of matching single or double quotes, like C++.
func Unquote(s string) string {
	s = Trim(s)
	if len(s) >= 2 && ((s[0] == '\'' && s[len(s)-1] == '\'') ||
		(s[0] == '"' && s[len(s)-1] == '"')) {
		return s[1 : len(s)-1]
	}
	return s
}

// StripComment cuts a line at the first '!' outside quotes, like C++.
func StripComment(line string) string {
	var single, dbl bool
	for i := 0; i < len(line); i++ {
		switch {
		case line[i] == '\'' && !dbl:
			single = !single
		case line[i] == '"' && !single:
			dbl = !dbl
		case line[i] == '!' && !single && !dbl:
			return line[:i]
		}
	}
	return line
}

// stripKeyValueComment cuts at '!' or '#' outside quotes, like C++.
func stripKeyValueComment(line string) string {
	var single, dbl bool
	for i := 0; i < len(line); i++ {
		switch {
		case line[i] == '\'' && !dbl:
			single = !single
		case line[i] == '"' && !single:
			dbl = !dbl
		case (line[i] == '!' || line[i] == '#') && !single && !dbl:
			return line[:i]
		}
	}
	return line
}

// Number parses a Fortran-style number (D/d exponents) and rejects
// non-finite values, like C++ number().
func Number(s string) (float64, error) {
	s = Unquote(Trim(s))
	s = strings.ReplaceAll(s, "D", "E")
	s = strings.ReplaceAll(s, "d", "e")
	x, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, fmt.Errorf("Bad number: %s", s)
	}
	if math.IsInf(x, 0) || math.IsNaN(x) {
		return 0, fmt.Errorf("Non-finite number: %s", s)
	}
	return x, nil
}

// ReadText reads a whole file as a string, like C++ read_text().
func ReadText(p string) (string, error) {
	data, err := os.ReadFile(p)
	if err != nil {
		return "", fmt.Errorf("Cannot read %s", p)
	}
	return string(data), nil
}

// WriteText writes s to p, creating parent directories, like C++ write_text().
func WriteText(p, s string) error {
	if dir := filepath.Dir(p); dir != "" && dir != "." {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			return fmt.Errorf("Cannot write %s", p)
		}
	}
	if err := os.WriteFile(p, []byte(s), 0o644); err != nil {
		return fmt.Errorf("Cannot write %s", p)
	}
	return nil
}

// DisplayPath renders p relative to the current working directory when
// possible, otherwise as an absolute normalized path, like C++ display_path().
func DisplayPath(p string) string {
	cwd, err := os.Getwd()
	if err != nil {
		return filepath.Clean(p)
	}
	var abs string
	if filepath.IsAbs(p) {
		abs = p
	} else {
		abs = filepath.Join(cwd, p)
	}
	abs = filepath.Clean(abs)
	rel, err := filepath.Rel(cwd, abs)
	if err == nil && rel != "" {
		return filepath.Clean(rel)
	}
	return abs
}

// nonFinite renders inf/nan like libstdc++ iostream: lowercase "inf",
// "-inf", "nan" (uppercase exponent variants handled by the callers).
func nonFinite(x float64) (string, bool) {
	switch {
	case math.IsNaN(x):
		return "nan", true
	case math.IsInf(x, 1):
		return "inf", true
	case math.IsInf(x, -1):
		return "-inf", true
	}
	return "", false
}

// FormatGeneral reproduces C++ iostream defaultfloat output at the given
// precision. libstdc++ implements defaultfloat as printf %g (trailing zeros
// removed, exponent threshold < -4 or >= precision), which is exactly Go's
// 'g' format; both round correctly to the same decimal digits.
func FormatGeneral(x float64, prec int) string {
	if s, ok := nonFinite(x); ok {
		return s
	}
	return strconv.FormatFloat(x, 'g', prec, 64)
}

// FormatSci reproduces C++ scientific output at the given precision, which
// counts digits after the decimal point (e.g. setprecision(15) -> 15
// decimals), like printf %e.
func FormatSci(x float64, prec int) string {
	if s, ok := nonFinite(x); ok {
		return s
	}
	return strconv.FormatFloat(x, 'e', prec, 64)
}

// FormatSciUpper is FormatSci with an uppercase exponent (std::uppercase).
func FormatSciUpper(x float64, prec int) string {
	if s, ok := nonFinite(x); ok {
		return strings.ToUpper(s)
	}
	return strconv.FormatFloat(x, 'E', prec, 64)
}

// FormatFixed reproduces C++ fixed output with prec decimal places.
func FormatFixed(x float64, prec int) string {
	if s, ok := nonFinite(x); ok {
		return s
	}
	return strconv.FormatFloat(x, 'f', prec, 64)
}

// Right pads s to width with spaces (C++ setw right-alignment default).
func Right(s string, width int) string {
	if len(s) >= width {
		return s
	}
	return strings.Repeat(" ", width-len(s)) + s
}

// Left pads s to width with trailing spaces (C++ setw with std::left).
func Left(s string, width int) string {
	if len(s) >= width {
		return s
	}
	return s + strings.Repeat(" ", width-len(s))
}

// Vec3 is a fixed-size three-component vector, like std::array<double, 3>.
type Vec3 [3]float64

// Config is a lowercased key=value map parsed from a plain text file.
type Config struct {
	v map[string]string
}

func loadConfigImpl(p string, allowNamelist bool) (*Config, error) {
	text, err := ReadText(p)
	if err != nil {
		return nil, err
	}
	c := &Config{v: make(map[string]string)}
	for _, raw := range strings.Split(text, "\n") {
		var line string
		if allowNamelist {
			line = Trim(StripComment(raw))
		} else {
			line = Trim(stripKeyValueComment(raw))
		}
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, "&") || line == "/" {
			if allowNamelist {
				continue
			}
			return nil, fmt.Errorf("Namelist syntax is not accepted in %s: %s", p, line)
		}
		eq := strings.IndexByte(line, '=')
		if eq < 0 {
			return nil, fmt.Errorf("Malformed assignment in %s: %s", p, line)
		}
		key := Lower(Trim(line[:eq]))
		val := Trim(line[eq+1:])
		if len(val) > 0 && val[len(val)-1] == ',' {
			val = val[:len(val)-1]
		}
		if key == "" || val == "" {
			return nil, fmt.Errorf("Malformed assignment in %s: %s", p, line)
		}
		if _, dup := c.v[key]; dup {
			return nil, fmt.Errorf("Duplicate parameter in %s: %s", p, key)
		}
		c.v[key] = Trim(val)
	}
	return c, nil
}

// Load parses a plain key=value file, like Config::load.
func Load(p string) (*Config, error) {
	return loadConfigImpl(p, false)
}

// LoadQENamelist parses a QE namelist-style file (&name ... /), like
// Config::load_qe_namelist.
func LoadQENamelist(p string) (*Config, error) {
	return loadConfigImpl(p, true)
}

// Has reports whether key is present (case-insensitive).
func (c *Config) Has(key string) bool {
	_, ok := c.v[Lower(key)]
	return ok
}

// Get returns the unquoted value for key, or d when absent.
func (c *Config) Get(key, d string) string {
	if v, ok := c.v[Lower(key)]; ok {
		return Unquote(v)
	}
	return d
}

// Real returns the numeric value for key, or d when absent.
func (c *Config) Real(key string, d float64) (float64, error) {
	if !c.Has(key) {
		return d, nil
	}
	return Number(c.v[Lower(key)])
}

// Integer returns the integral value for key, or d when absent; non-integral
// values are rejected, like Config::integer.
func (c *Config) Integer(key string, d int) (int, error) {
	if !c.Has(key) {
		return d, nil
	}
	x, err := Number(c.v[Lower(key)])
	if err != nil {
		return 0, err
	}
	rounded := math.Round(x)
	if math.Abs(x-rounded) > 1.0e-10 {
		return 0, fmt.Errorf("Expected integer for %s: %s", key, c.v[Lower(key)])
	}
	return int(rounded), nil
}

// RequireOnly rejects any key outside the allowed set, like
// Config::require_only.
func (c *Config) RequireOnly(allowed map[string]bool, context string) error {
	keys := make([]string, 0, len(c.v))
	for key := range c.v {
		keys = append(keys, key)
	}
	sort.Strings(keys)
	for _, key := range keys {
		if !allowed[key] {
			return fmt.Errorf("Unknown parameter in %s: %s", context, key)
		}
	}
	return nil
}

// LogicalValue interprets a Fortran-style logical, like logical_value().
func LogicalValue(c *Config, key string, required bool) (bool, error) {
	if !c.Has(key) {
		if required {
			return false, fmt.Errorf("Missing required parameter in fdvib.in: %s", key)
		}
		return false, nil
	}
	switch value := Lower(c.Get(key, "")); value {
	case ".true.", "true":
		return true, nil
	case ".false.", "false":
		return false, nil
	default:
		return false, fmt.Errorf("Expected logical value for %s: %s", key, c.Get(key, ""))
	}
}

// IntegerList parses a comma/semicolon/space separated list of integers,
// like integer_list().
func IntegerList(s string) ([]int, error) {
	original := s
	s = strings.Map(func(r rune) rune {
		if r == ',' || r == ';' {
			return ' '
		}
		return r
	}, s)
	fields := strings.Fields(s)
	out := make([]int, 0, len(fields))
	for _, f := range fields {
		x, err := strconv.Atoi(f)
		if err != nil {
			return nil, fmt.Errorf("Bad integer list: %s", original)
		}
		out = append(out, x)
	}
	return out, nil
}
