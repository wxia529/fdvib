// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

#include "fdvib.hpp"

#include <array>
#include <stdexcept>

namespace fdvib {

namespace {

const std::array<const char *, 119> &element_symbols() {
    static const std::array<const char *, 119> symbols = {
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
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    };
    return symbols;
}

int atomic_number(const std::string &symbol) {
    const auto &symbols = element_symbols();
    for (std::size_t z = 1; z < symbols.size(); ++z)
        if (symbol == symbols[z]) return static_cast<int>(z);
    return 0;
}

} // namespace

bool is_standard_element(const std::string &symbol) {
    return atomic_number(symbol) != 0;
}

std::string standard_element_symbol(const std::string &label, const std::string &context) {
    if (is_standard_element(label)) return label;
    if (label.size() >= 2) {
        const auto two = label.substr(0, 2);
        if (is_standard_element(two)) return two;
    }
    if (!label.empty()) {
        const auto one = label.substr(0, 1);
        if (is_standard_element(one)) return one;
    }
    throw std::runtime_error("Cannot map species label to a standard element for " + context + ": " + label);
}

int atomic_number_from_label(const std::string &label, const std::string &context) {
    return atomic_number(standard_element_symbol(label, context));
}

} // namespace fdvib
