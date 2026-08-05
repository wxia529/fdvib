// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

// Package units holds the physical constants used throughout FDVIB.
// Values mirror fdvib.hpp (CODATA 2022). AMU_RY = 1 / (2 m_e[u]) for QE
// Rydberg atomic units.
package units

const (
	BohrToAng = 0.529177210544
	AmuRy     = 911.4442431390707
	CMToEV    = 1.2398419843320026e-4
	KBEV      = 8.617333262145179e-5
	KBSI      = 1.380649e-23
	HSI       = 6.62607015e-34
	AmuKg     = 1.66053906892e-27
	AtmPa     = 101325.0
	PI        = 3.141592653589793238462643383279502884
)
