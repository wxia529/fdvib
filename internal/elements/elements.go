// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package elements provides the periodic-table symbol table and species-label
// mapping, mirroring elements.cpp.
package elements

import "fmt"

// elementSymbols are the standard element symbols indexed by atomic number
// (1-based; index 0 is unused), matching elements.cpp.
var elementSymbols = [119]string{
	"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
	"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
	"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
	"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
	"Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
	"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
	"Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
	"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
	"Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
	"Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
	"Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}

// atomicNumber returns the atomic number of a standard element symbol,
// or 0 when the symbol is not standard.
func atomicNumber(symbol string) int {
	for z := 1; z < len(elementSymbols); z++ {
		if symbol == elementSymbols[z] {
			return z
		}
	}
	return 0
}

// IsStandardElement reports whether symbol is a standard element symbol.
func IsStandardElement(symbol string) bool {
	return atomicNumber(symbol) != 0
}

// StandardElementSymbol reduces a species label (e.g. "C1", "O_ads") to the
// leading standard element symbol, like standard_element_symbol().
func StandardElementSymbol(label, context string) (string, error) {
	if IsStandardElement(label) {
		return label, nil
	}
	if len(label) >= 2 {
		two := label[:2]
		if IsStandardElement(two) {
			return two, nil
		}
	}
	if len(label) > 0 {
		one := label[:1]
		if IsStandardElement(one) {
			return one, nil
		}
	}
	return "", fmt.Errorf("Cannot map species label to a standard element for %s: %s", context, label)
}

// AtomicNumberFromLabel returns the atomic number of the standard element
// mapped from label, like atomic_number_from_label().
func AtomicNumberFromLabel(label, context string) (int, error) {
	symbol, err := StandardElementSymbol(label, context)
	if err != nil {
		return 0, err
	}
	return atomicNumber(symbol), nil
}
