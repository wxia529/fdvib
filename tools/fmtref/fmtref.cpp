// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia
//
// Reference generator for the C++ iostream formatting helpers in
// internal/config (FormatGeneral/FormatSci/FormatSciUpper/FormatFixed).
//
// Regenerate the byte-exact reference after any intentional change to the
// formatting semantics:
//
//	g++ -O2 -o fmtref fmtref.cpp
//	./fmtref > ../../internal/config/testdata/fmt_cpp.txt
//
// The Go helpers must keep reproducing this output byte-for-byte
// (TestFormatMatchesCPP in internal/config). Change the generator first,
// regenerate the reference, then adjust the Go helpers.

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

int main() {
    std::vector<double> vals = {
        0.0, -0.0, 1.0, -1.0, 1.5, 0.5, 0.001, 123.45, 999.9999999999999,
        1e5, 1e-5, 1.23e-7, -3.5, 123456789.0, 0.1, 0.01, 1.0/3.0, 2.5e10,
        3e-10, 1e15, 123456789012345.0, 0.000123456789, 6.02214076e23,
        5.0e-324, 2.2250738585072014e-308, 1.7976931348623157e308, 1000.0,
        1000000000000000.0, 123456789.123456789, 0.0001, 0.00009999999999999999
    };
    const int precs[] = {6, 7, 9, 10, 12, 17};
    for (double x : vals) {
        for (int p : precs) {
            std::ostringstream o; o << std::setprecision(p) << x;
            std::cout << "G" << p << ":" << o.str() << "\n";
        }
        for (int p : precs) {
            std::ostringstream o; o << std::scientific << std::setprecision(p) << x;
            std::cout << "S" << p << ":" << o.str() << "\n";
        }
        for (int p : precs) {
            std::ostringstream o; o << std::fixed << std::setprecision(p) << x;
            std::cout << "F" << p << ":" << o.str() << "\n";
        }
    }
}
