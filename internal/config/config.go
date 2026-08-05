// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package config provides text utilities, the strict key=value configuration
// parser, and the numeric formatting helpers for the stable output formats
// consumed by external tools (QE's dynmat.x, Shermo, Molden viewers) and by
// FDVIB's own restart markers. The exact byte layout is a frozen contract.
//
// Formatting semantics that differ from strconv/fmt:
//
//   - FormatSci and FormatFixed treat prec as digits after the decimal point
//     (std::setprecision under scientific/fixed), not significant digits.
//   - FormatGeneral is printf %g (std::defaultfloat): trailing zeros are
//     stripped and exponential form is used when the exponent is < -4 or
//     >= prec.
//   - Non-finite values render as "inf"/"-inf"/"nan" (uppercase with
//     FormatSciUpper), matching libstdc++.
//
// Do not replace these helpers with direct strconv calls. The byte-exact
// baseline testdata/fmt_cpp.txt is frozen; changing the output format is a
// contract change for external consumers and requires updating the baseline
// and TestFormatMatchesCPP together.
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

// Unquote removes one layer of matching single or double quotes, like C++.
func Unquote(s string) string {
	s = strings.Trim(s, " \t\r\n")
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

// ParseNumber parses a Fortran-style number (D/d exponents) and rejects
// non-finite values, like C++ number().
func ParseNumber(s string) (float64, error) {
	s = Unquote(strings.Trim(s, " \t\r\n"))
	s = strings.ReplaceAll(s, "D", "E")
	s = strings.ReplaceAll(s, "d", "e")
	x, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, fmt.Errorf("invalid number %s", s)
	}
	if math.IsInf(x, 0) || math.IsNaN(x) {
		return 0, fmt.Errorf("non-finite number %s", s)
	}
	return x, nil
}

// ReadText reads a whole file as a string, like C++ read_text().
func ReadText(p string) (string, error) {
	data, err := os.ReadFile(p)
	if err != nil {
		return "", fmt.Errorf("cannot read %s", p)
	}
	return string(data), nil
}

// WriteText writes s to p, creating parent directories, like C++ write_text().
func WriteText(p, s string) error {
	if dir := filepath.Dir(p); dir != "" && dir != "." {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			return fmt.Errorf("cannot write %s", p)
		}
	}
	if err := os.WriteFile(p, []byte(s), 0o644); err != nil {
		return fmt.Errorf("cannot write %s", p)
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

// FixedField renders v with prec decimal places (C++ fixed + setprecision),
// right-padded to width (C++ setw). It corresponds to
// `setw(width) << fixed << setprecision(prec) << v`.
func FixedField(v float64, prec, width int) string {
	return fmt.Sprintf("%*s", width, FormatFixed(v, prec))
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
			line = strings.Trim(StripComment(raw), " \t\r\n")
		} else {
			line = strings.Trim(stripKeyValueComment(raw), " \t\r\n")
		}
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, "&") || line == "/" {
			if allowNamelist {
				continue
			}
			return nil, fmt.Errorf("namelist syntax is not accepted in %s: %s", p, line)
		}
		eq := strings.IndexByte(line, '=')
		if eq < 0 {
			return nil, fmt.Errorf("malformed assignment in %s: %s", p, line)
		}
		key := strings.ToLower(strings.Trim(line[:eq], " \t\r\n"))
		val := strings.Trim(line[eq+1:], " \t\r\n")
		if len(val) > 0 && val[len(val)-1] == ',' {
			val = val[:len(val)-1]
		}
		if key == "" || val == "" {
			return nil, fmt.Errorf("malformed assignment in %s: %s", p, line)
		}
		if _, dup := c.v[key]; dup {
			return nil, fmt.Errorf("duplicate parameter in %s: %s", p, key)
		}
		c.v[key] = strings.Trim(val, " \t\r\n")
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
	_, ok := c.v[strings.ToLower(key)]
	return ok
}

// Get returns the unquoted value for key, or d when absent.
func (c *Config) Get(key, d string) string {
	if v, ok := c.v[strings.ToLower(key)]; ok {
		return Unquote(v)
	}
	return d
}

// Real returns the numeric value for key, or d when absent.
func (c *Config) Real(key string, d float64) (float64, error) {
	if !c.Has(key) {
		return d, nil
	}
	return ParseNumber(c.v[strings.ToLower(key)])
}

// Integer returns the integral value for key, or d when absent; non-integral
// values are rejected, like Config::integer.
func (c *Config) Integer(key string, d int) (int, error) {
	if !c.Has(key) {
		return d, nil
	}
	x, err := ParseNumber(c.v[strings.ToLower(key)])
	if err != nil {
		return 0, err
	}
	rounded := math.Round(x)
	if math.Abs(x-rounded) > 1.0e-10 {
		return 0, fmt.Errorf("expected integer for %s: %s", key, c.v[strings.ToLower(key)])
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
			return fmt.Errorf("unknown parameter in %s: %s", context, key)
		}
	}
	return nil
}

// ParseLogical interprets a Fortran-style logical, like logical_value().
func ParseLogical(c *Config, key string, required bool) (bool, error) {
	if !c.Has(key) {
		if required {
			return false, fmt.Errorf("missing required parameter in fdvib.in: %s", key)
		}
		return false, nil
	}
	switch value := strings.ToLower(c.Get(key, "")); value {
	case ".true.", "true":
		return true, nil
	case ".false.", "false":
		return false, nil
	default:
		return false, fmt.Errorf("expected logical value for %s: %s", key, c.Get(key, ""))
	}
}

// ParseIntList parses a comma/semicolon/space separated list of integers,
// like integer_list().
func ParseIntList(s string) ([]int, error) {
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
			return nil, fmt.Errorf("invalid integer list %s", original)
		}
		out = append(out, x)
	}
	return out, nil
}

// SplitLines splits text into lines with std::getline semantics: a trailing
// newline does not produce an extra empty line.
func SplitLines(text string) []string {
	if text == "" {
		return nil
	}
	lines := strings.Split(text, "\n")
	if lines[len(lines)-1] == "" {
		lines = lines[:len(lines)-1]
	}
	return lines
}

// IstreamDouble mimics C++ istream >> double, which parses the longest
// strtod-compatible prefix and leaves the rest unconsumed (strtod stops at
// Fortran D/d exponents, so "0.757D0" parses as 0.757 with remainder "D0").
// The returned ok is false when no numeric prefix exists.
func IstreamDouble(s string) (float64, string, bool) {
	s = strings.TrimLeft(s, " \t")
	for i := len(s); i >= 1; i-- {
		if x, err := strconv.ParseFloat(s[:i], 64); err == nil {
			return x, s[i:], true
		}
	}
	return 0, s, false
}
