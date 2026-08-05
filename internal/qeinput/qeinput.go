// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package qeinput parses Quantum ESPRESSO pw.x input files and rewrites them
// for reference and displaced attempts, mirroring qe_input.cpp.
package qeinput

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/wxia529/fdvib/internal/config"
	"github.com/wxia529/fdvib/internal/units"
)

// Species is an ATOMIC_SPECIES entry, like struct Species.
type Species struct {
	Symbol string
	Pseudo string
	Mass   float64
}

// Atom is an ATOMIC_POSITIONS entry, like struct Atom.
type Atom struct {
	Symbol string
	// R is the Cartesian position in angstrom.
	R config.Vec3
	// InputR keeps the coordinates exactly as written on the card.
	InputR config.Vec3
	Extra  []string
	// Type is the 1-based species index.
	Type int
}

// QEInput is the parsed and normalized scf.in, like struct QEInput.
type QEInput struct {
	Text      string
	CleanText string
	Prefix    string
	// Lines preserves the original card text, each line ending in "\n".
	Lines         []string
	Nat, Ntyp     int
	PosHeader     int
	PosStart      int
	CellHeader    int
	PositionsUnit string
	CellUnit      string
	AlatAngstrom  float64
	Species       []Species
	Atoms         []Atom
	// Cell is the 3x3 cell in angstrom.
	Cell [3]config.Vec3
}

var (
	natRe           = regexp.MustCompile(`(?i)\bnat\s*=\s*(\d+)`)
	ntypRe          = regexp.MustCompile(`(?i)\bntyp\s*=\s*(\d+)`)
	ibravRe         = regexp.MustCompile(`(?i)\bibrav\s*=\s*([-+]?\d+)`)
	scfRe           = regexp.MustCompile(`(?i)\bcalculation\s*=\s*['"]scf['"]`)
	tprnforRe       = regexp.MustCompile(`(?i)\btprnfor\s*=\s*\.true\.`)
	startingpotRe   = regexp.MustCompile(`(?i)\bstartingpot\s*=\s*['"]([^'"]+)['"]`)
	prefixRe        = regexp.MustCompile(`(?i)\bprefix\s*=\s*['"]([^'"]+)['"]`)
	speciesCardRe   = regexp.MustCompile(`(?i)^ATOMIC_SPECIES\b`)
	positionsCardRe = regexp.MustCompile(`(?i)^ATOMIC_POSITIONS\b`)
	cellCardRe      = regexp.MustCompile(`(?i)^CELL_PARAMETERS\b`)
)

func cardUnit(line, card string) string {
	re := regexp.MustCompile(`(?i)^\s*` + card + `\b\s*(?:[({]\s*([A-Za-z_]+)\s*[)}]|([A-Za-z_]+))?`)
	m := re.FindStringSubmatch(line)
	if m == nil {
		return ""
	}
	unit := m[1]
	if unit == "" {
		unit = m[2]
	}
	return config.Lower(unit)
}

func scalarParameter(text, name string) (float64, bool) {
	re := regexp.MustCompile(`(?i)\b` + name + `\s*=\s*([-+0-9.EeDd]+)`)
	m := re.FindStringSubmatch(text)
	if m == nil {
		return 0.0, false
	}
	x, err := config.Number(m[1])
	if err != nil {
		return 0.0, false
	}
	return x, true
}

func alatAngstrom(text string) float64 {
	if a, ok := scalarParameter(text, "A"); ok && a > 0.0 {
		return a
	}
	if celldm1, ok := scalarParameter(text, `celldm\s*\(\s*1\s*\)`); ok && celldm1 > 0.0 {
		return celldm1 * units.BohrToAng
	}
	return 0.0
}

func scaled(v config.Vec3, factor float64) config.Vec3 {
	return config.Vec3{v[0] * factor, v[1] * factor, v[2] * factor}
}

func crystalToCart(frac config.Vec3, cell [3]config.Vec3) config.Vec3 {
	var out config.Vec3
	for i := 0; i < 3; i++ {
		for k := 0; k < 3; k++ {
			out[k] += frac[i] * cell[i][k]
		}
	}
	return out
}

func cartToCrystalDelta(cell [3]config.Vec3, cart config.Vec3) (config.Vec3, error) {
	a00, a01, a02 := cell[0][0], cell[1][0], cell[2][0]
	a10, a11, a12 := cell[0][1], cell[1][1], cell[2][1]
	a20, a21, a22 := cell[0][2], cell[1][2], cell[2][2]
	det := a00*(a11*a22-a12*a21) - a01*(a10*a22-a12*a20) +
		a02*(a10*a21-a11*a20)
	if math.Abs(det) < 1.0e-14 {
		return config.Vec3{}, fmt.Errorf("CELL_PARAMETERS matrix is singular")
	}
	var out config.Vec3
	out[0] = ((a11*a22-a12*a21)*cart[0] + (a02*a21-a01*a22)*cart[1] +
		(a01*a12-a02*a11)*cart[2]) / det
	out[1] = ((a12*a20-a10*a22)*cart[0] + (a00*a22-a02*a20)*cart[1] +
		(a02*a10-a00*a12)*cart[2]) / det
	out[2] = ((a10*a21-a11*a20)*cart[0] + (a01*a20-a00*a21)*cart[1] +
		(a00*a11-a01*a10)*cart[2]) / det
	return out, nil
}

func supportedPositionUnit(unit string) bool {
	return unit == "angstrom" || unit == "bohr" || unit == "alat" || unit == "crystal"
}

func supportedCellUnit(unit string) bool {
	return unit == "angstrom" || unit == "bohr" || unit == "alat"
}

// ParseQeInput reads and validates a pw.x input file, like parse_qe_input().
func ParseQeInput(p string) (*QEInput, error) {
	text, err := config.ReadText(p)
	if err != nil {
		return nil, err
	}
	q := &QEInput{Text: text, Prefix: "pwscf",
		PositionsUnit: "angstrom", CellUnit: "angstrom",
		PosHeader: -1, PosStart: -1, CellHeader: -1}
	for _, line := range config.SplitLines(text) {
		q.Lines = append(q.Lines, line+"\n")
		q.CleanText += config.StripComment(line) + "\n"
	}

	mnat := natRe.FindStringSubmatch(q.CleanText)
	if mnat == nil {
		return nil, fmt.Errorf("cannot find nat in %s", p)
	}
	mntyp := ntypRe.FindStringSubmatch(q.CleanText)
	if mntyp == nil {
		return nil, fmt.Errorf("cannot find ntyp in %s", p)
	}
	q.Nat, _ = strconv.Atoi(mnat[1])
	q.Ntyp, _ = strconv.Atoi(mntyp[1])
	if q.Nat <= 0 || q.Ntyp <= 0 {
		return nil, fmt.Errorf("nat and ntyp must be positive")
	}
	mib := ibravRe.FindStringSubmatch(q.CleanText)
	if mib == nil {
		return nil, fmt.Errorf("cannot find ibrav in %s", p)
	}
	ibrav, _ := strconv.Atoi(mib[1])
	if ibrav != 0 {
		return nil, fmt.Errorf("fdvib requires ibrav=0")
	}
	if !scfRe.MatchString(q.CleanText) {
		return nil, fmt.Errorf("scf.in must contain calculation='scf'")
	}
	if !tprnforRe.MatchString(q.CleanText) {
		return nil, fmt.Errorf("scf.in must contain tprnfor=.true.")
	}
	for _, m := range startingpotRe.FindAllStringSubmatch(q.CleanText, -1) {
		if config.Lower(m[1]) == "file" {
			return nil, fmt.Errorf("scf.in must not set startingpot='file'; fdvib manages the reference density")
		}
	}
	if m := prefixRe.FindStringSubmatch(q.CleanText); m != nil {
		q.Prefix = m[1]
	}
	if q.Prefix == "" || filepath.Base(q.Prefix) != q.Prefix ||
		q.Prefix == "." || q.Prefix == ".." {
		return nil, fmt.Errorf("qe prefix must be a non-empty filename prefix without directories")
	}

	q.AlatAngstrom = alatAngstrom(q.CleanText)
	sp := -1
	for i, raw := range q.Lines {
		s := config.Trim(raw)
		switch {
		case speciesCardRe.MatchString(s):
			sp = i
		case positionsCardRe.MatchString(s):
			q.PositionsUnit = cardUnit(s, "ATOMIC_POSITIONS")
			if q.PositionsUnit == "" {
				return nil, fmt.Errorf("ATOMIC_POSITIONS must specify units: angstrom, bohr, alat, or crystal")
			}
			if q.PositionsUnit == "crystal_sg" {
				return nil, fmt.Errorf("ATOMIC_POSITIONS crystal_sg is not supported; use explicit crystal coordinates")
			}
			if !supportedPositionUnit(q.PositionsUnit) {
				return nil, fmt.Errorf("unsupported ATOMIC_POSITIONS unit %s", q.PositionsUnit)
			}
			q.PosHeader = i
			q.PosStart = i + 1
		case cellCardRe.MatchString(s):
			q.CellUnit = cardUnit(s, "CELL_PARAMETERS")
			if q.CellUnit == "" {
				return nil, fmt.Errorf("CELL_PARAMETERS must specify units: angstrom, bohr, or alat")
			}
			if !supportedCellUnit(q.CellUnit) {
				return nil, fmt.Errorf("unsupported CELL_PARAMETERS unit %s", q.CellUnit)
			}
			q.CellHeader = i
		}
	}
	if sp < 0 || q.PosStart < 0 || q.CellHeader < 0 {
		return nil, fmt.Errorf("require ATOMIC_SPECIES, ATOMIC_POSITIONS, and CELL_PARAMETERS in scf.in")
	}
	lineCount := len(q.Lines)
	if sp+q.Ntyp >= lineCount {
		return nil, fmt.Errorf("ATOMIC_SPECIES block is shorter than ntyp")
	}
	if q.PosStart+q.Nat > lineCount {
		return nil, fmt.Errorf("ATOMIC_POSITIONS block is shorter than nat")
	}
	if q.CellHeader+3 >= lineCount {
		return nil, fmt.Errorf("CELL_PARAMETERS block is incomplete")
	}
	if (q.PositionsUnit == "alat" || q.CellUnit == "alat") && q.AlatAngstrom <= 0.0 {
		return nil, fmt.Errorf("alat units require A or celldm(1) in &SYSTEM")
	}

	for i := 0; i < q.Ntyp; i++ {
		fields := strings.Fields(q.Lines[sp+1+i])
		if len(fields) < 3 {
			return nil, fmt.Errorf("malformed ATOMIC_SPECIES line")
		}
		var x Species
		x.Symbol = fields[0]
		mass, _, ok := config.IstreamDouble(fields[1])
		if !ok {
			return nil, fmt.Errorf("malformed ATOMIC_SPECIES line")
		}
		x.Mass = mass
		x.Pseudo = fields[2]
		if !(x.Mass > 0.0) || math.IsInf(x.Mass, 0) || math.IsNaN(x.Mass) {
			return nil, fmt.Errorf("atomic mass must be positive")
		}
		q.Species = append(q.Species, x)
	}
	for i := 0; i < q.Nat; i++ {
		fields := strings.Fields(q.Lines[q.PosStart+i])
		if len(fields) < 4 {
			return nil, fmt.Errorf("malformed ATOMIC_POSITIONS line")
		}
		var a Atom
		a.Symbol = fields[0]
		for k := 0; k < 3; k++ {
			// A non-empty remainder means the next istream read would have
			// failed, so the whole line is rejected like C++.
			x, rest, ok := config.IstreamDouble(fields[1+k])
			if !ok || rest != "" {
				return nil, fmt.Errorf("malformed ATOMIC_POSITIONS line")
			}
			a.InputR[k] = x
		}
		a.Extra = append(a.Extra, fields[4:]...)
		typeIndex := -1
		for si, s := range q.Species {
			if s.Symbol == a.Symbol {
				typeIndex = si
				break
			}
		}
		if typeIndex < 0 {
			return nil, fmt.Errorf("unknown species %s", a.Symbol)
		}
		a.Type = typeIndex + 1
		q.Atoms = append(q.Atoms, a)
	}
	for i := 0; i < 3; i++ {
		fields := strings.Fields(q.Lines[q.CellHeader+1+i])
		if len(fields) < 3 {
			return nil, fmt.Errorf("malformed CELL_PARAMETERS")
		}
		for k := 0; k < 3; k++ {
			x, rest, ok := config.IstreamDouble(fields[k])
			if !ok || rest != "" {
				return nil, fmt.Errorf("malformed CELL_PARAMETERS")
			}
			q.Cell[i][k] = x
		}
	}
	var cellFactor float64
	switch q.CellUnit {
	case "angstrom":
		cellFactor = 1.0
	case "bohr":
		cellFactor = units.BohrToAng
	default: // alat
		cellFactor = q.AlatAngstrom
	}
	for i := range q.Cell {
		q.Cell[i] = scaled(q.Cell[i], cellFactor)
	}
	for i := range q.Atoms {
		switch q.PositionsUnit {
		case "angstrom":
			q.Atoms[i].R = q.Atoms[i].InputR
		case "bohr":
			q.Atoms[i].R = scaled(q.Atoms[i].InputR, units.BohrToAng)
		case "alat":
			q.Atoms[i].R = scaled(q.Atoms[i].InputR, q.AlatAngstrom)
		case "crystal":
			q.Atoms[i].R = crystalToCart(q.Atoms[i].InputR, q.Cell)
		}
	}
	return q, nil
}

// FormatPosition renders one atom card line with C++ setw/setprecision
// semantics, like format_position().
func FormatPosition(a *Atom, coord config.Vec3) string {
	var b strings.Builder
	b.WriteString(config.Left(a.Symbol, 4))
	for _, v := range coord {
		b.WriteString(" ")
		b.WriteString(config.Right(config.FormatFixed(v, 10), 18))
	}
	for _, x := range a.Extra {
		b.WriteString("  ")
		b.WriteString(x)
	}
	b.WriteString("\n")
	return b.String()
}

var (
	outdirRe          = regexp.MustCompile(`(?i)^\s*outdir\s*=`)
	wfcdirRe          = regexp.MustCompile(`(?i)^\s*wfcdir\s*=`)
	pseudoRe          = regexp.MustCompile(`(?i)^\s*pseudo_dir\s*=\s*(['"])([^'"]+)['"]`)
	diskIORe          = regexp.MustCompile(`(?i)^\s*disk_io\s*=`)
	startingpotLineRe = regexp.MustCompile(`(?i)^\s*startingpot\s*=`)
	controlRe         = regexp.MustCompile(`(?i)^&CONTROL\b`)
	electronsRe       = regexp.MustCompile(`(?i)^&ELECTRONS\b`)
)

// AttemptLines rewrites outdir/wfcdir/pseudo_dir for an isolated attempt
// directory, like attempt_lines().
func AttemptLines(q *QEInput, outdir string, sourceDir, runDir string) ([]string, error) {
	lines := make([]string, len(q.Lines))
	copy(lines, q.Lines)
	found := false
	for i, line := range lines {
		indent := line[:len(line)-len(strings.TrimLeft(line, " \t"))]
		switch {
		case outdirRe.MatchString(line):
			lines[i] = indent + "outdir = '" + outdir + "',\n"
			found = true
		case wfcdirRe.MatchString(line):
			lines[i] = indent + "wfcdir = '" + outdir + "',\n"
		default:
			m := pseudoRe.FindStringSubmatch(line)
			if m == nil {
				continue
			}
			original := m[2]
			value := original
			homeRelative := original == "~" || strings.HasPrefix(original, "~/") ||
				original == "$HOME" || strings.HasPrefix(original, "$HOME/") ||
				original == "${HOME}" || strings.HasPrefix(original, "${HOME}/")
			switch {
			case homeRelative:
				home := os.Getenv("HOME")
				if home == "" {
					return nil, fmt.Errorf("pseudo_dir uses the home directory but HOME is not set")
				}
				var suffix string
				switch {
				case strings.HasPrefix(original, "~/"):
					suffix = original[2:]
				case strings.HasPrefix(original, "$HOME/"):
					suffix = original[6:]
				case strings.HasPrefix(original, "${HOME}/"):
					suffix = original[8:]
				}
				if filepath.IsAbs(suffix) {
					value = suffix
				} else {
					value = filepath.Join(home, suffix)
				}
			case !filepath.IsAbs(value) && original != "" &&
				original[0] != '~' && original[0] != '$':
				target := filepath.Clean(filepath.Join(sourceDir, value))
				rel, err := filepath.Rel(filepath.Clean(runDir), target)
				if err != nil || rel == "" {
					value = "."
				} else {
					value = rel
				}
			}
			lines[i] = indent + "pseudo_dir = '" + filepath.ToSlash(value) + "',\n"
		}
	}
	if !found {
		return nil, fmt.Errorf("cannot find outdir in scf.in")
	}
	return lines, nil
}

// ReferenceInput rewrites the scf.in for the reference attempt, like
// reference_input().
func ReferenceInput(q *QEInput, outdir string, sourceDir, runDir string) (string, error) {
	lines, err := AttemptLines(q, outdir, sourceDir, runDir)
	if err != nil {
		return "", err
	}
	return strings.Join(lines, ""), nil
}

// DisplacedInput rewrites the scf.in for a displaced attempt, like
// displaced_input().
func DisplacedInput(q *QEInput, atom, axis int, shift float64,
	outdir string, sourceDir, runDir string) (string, error) {
	lines, err := AttemptLines(q, outdir, sourceDir, runDir)
	if err != nil {
		return "", err
	}
	a := q.Atoms[atom]
	coord := a.InputR
	switch q.PositionsUnit {
	case "angstrom":
		coord[axis] += shift
	case "bohr":
		coord[axis] += shift / units.BohrToAng
	case "alat":
		coord[axis] += shift / q.AlatAngstrom
	case "crystal":
		var cart config.Vec3
		cart[axis] = shift
		delta, err := cartToCrystalDelta(q.Cell, cart)
		if err != nil {
			return "", err
		}
		for k := 0; k < 3; k++ {
			coord[k] += delta[k]
		}
	}
	lines[q.PosStart+atom] = FormatPosition(&a, coord)

	diskIOFound := false
	for i, line := range lines {
		if diskIORe.MatchString(line) {
			indent := line[:len(line)-len(strings.TrimLeft(line, " \t"))]
			lines[i] = indent + "disk_io = 'minimal',\n"
			diskIOFound = true
		}
	}
	if !diskIOFound {
		control := -1
		for i := range lines {
			if controlRe.MatchString(config.Trim(lines[i])) {
				control = i
				break
			}
		}
		if control < 0 {
			return "", fmt.Errorf("cannot find &CONTROL in scf.in")
		}
		end := control + 1
		for end < len(lines) && config.Trim(config.StripComment(lines[end])) != "/" {
			end++
		}
		if end == len(lines) {
			return "", fmt.Errorf("unterminated &CONTROL namelist")
		}
		lines = append(lines[:end], append([]string{"  disk_io = 'minimal',\n"}, lines[end:]...)...)
	}

	startingpotFound := false
	for i, l := range lines {
		if startingpotLineRe.MatchString(l) {
			indent := l[:len(l)-len(strings.TrimLeft(l, " \t"))]
			lines[i] = indent + "startingpot = 'file',\n"
			startingpotFound = true
		}
	}
	if !startingpotFound {
		electrons := -1
		for i := range lines {
			if electronsRe.MatchString(config.Trim(lines[i])) {
				electrons = i
				break
			}
		}
		if electrons >= 0 {
			end := electrons + 1
			for end < len(lines) && config.Trim(config.StripComment(lines[end])) != "/" {
				end++
			}
			if end == len(lines) {
				return "", fmt.Errorf("unterminated &ELECTRONS namelist")
			}
			lines = append(lines[:end], append([]string{"  startingpot = 'file',\n"}, lines[end:]...)...)
		} else {
			firstCard := q.PosHeader
			if q.CellHeader < firstCard {
				firstCard = q.CellHeader
			}
			for i := range lines {
				if speciesCardRe.MatchString(config.Trim(lines[i])) && i < firstCard {
					firstCard = i
				}
			}
			lines = append(lines[:firstCard],
				append([]string{"&ELECTRONS\n  startingpot = 'file',\n/\n"}, lines[firstCard:]...)...)
		}
	}
	return strings.Join(lines, ""), nil
}
