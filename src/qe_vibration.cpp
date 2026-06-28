// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

#include "fdvib.hpp"

#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace fdvib {

std::vector<Mode> read_qe_dynmat_modes(const fs::path &path, int nat) {
    if (nat <= 0)
        throw std::runtime_error("Invalid atom count while parsing modes: " + path.string());
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read " + path.string());
    const std::regex frequency(
        R"(freq\s*\(\s*\d+\s*\)\s*=\s*[-+0-9.EeDd]+\s*\[THz\]\s*=\s*([-+0-9.EeDd]+)\s*\[cm-1\])",
        std::regex::icase);
    const std::regex eigenvector(
        R"(\(\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*\))");
    std::vector<Mode> modes;
    std::string line;
    std::smatch match;
    while (std::getline(input, line)) {
        if (!std::regex_search(line, match, frequency)) continue;
        Mode mode;
        mode.freq = number(match[1]);
        mode.displacement.reserve(nat);
        for (int i = 0; i < nat; ++i) {
            if (!std::getline(input, line) || !std::regex_search(line, match, eigenvector))
                throw std::runtime_error("Bad eigenvector block in " + path.string());
            mode.displacement.push_back({number(match[1]), number(match[3]), number(match[5])});
        }
        modes.push_back(std::move(mode));
    }
    if (static_cast<int>(modes.size()) != 3 * nat)
        throw std::runtime_error("Expected " + std::to_string(3 * nat) + " modes, found " +
                                 std::to_string(modes.size()) + " in " + path.string());
    return modes;
}

DynGeometry read_qe_dyn_geometry(const fs::path &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read " + path.string());
    std::string line;
    for (int i = 0; i < 3; ++i)
        if (!std::getline(input, line))
            throw std::runtime_error("Incomplete dynG header: " + path.string());
    std::istringstream header(line);
    int ntyp, nat, ibrav;
    double alat;
    header >> ntyp >> nat >> ibrav >> alat;
    if (!header || ibrav != 0 || ntyp <= 0 || nat <= 0 || !(alat > 0.0))
        throw std::runtime_error("Unsupported dynG header");
    if (!std::getline(input, line))
        throw std::runtime_error("Incomplete dynG basis section: " + path.string());

    DynGeometry geometry;
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(input, line))
            throw std::runtime_error("Incomplete dynG basis vectors: " + path.string());
        std::istringstream basis(line);
        Vec3 vector{};
        if (!(basis >> vector[0] >> vector[1] >> vector[2]) ||
            !std::isfinite(vector[0]) || !std::isfinite(vector[1]) ||
            !std::isfinite(vector[2]))
            throw std::runtime_error("Bad basis vector in dynG: " + path.string());
        for (int k = 0; k < 3; ++k) geometry.cell_bohr[i][k] = vector[k] * alat;
    }
    const auto &a = geometry.cell_bohr[0];
    const auto &b = geometry.cell_bohr[1];
    const auto &c = geometry.cell_bohr[2];
    const double volume = a[0]*(b[1]*c[2]-b[2]*c[1])
                        - a[1]*(b[0]*c[2]-b[2]*c[0])
                        + a[2]*(b[0]*c[1]-b[1]*c[0]);
    if (!std::isfinite(volume) || std::abs(volume) < 1.0e-12)
        throw std::runtime_error("Singular basis vectors in dynG: " + path.string());

    std::vector<double> type_masses(ntyp);
    std::vector<std::string> type_symbols(ntyp);
    const std::regex quoted(R"('([^']*)')");
    for (int i = 0; i < ntyp; ++i) {
        if (!std::getline(input, line))
            throw std::runtime_error("Incomplete dynG species block: " + path.string());
        std::smatch match;
        if (!std::regex_search(line, match, quoted))
            throw std::runtime_error("Bad species symbol in dynG");
        type_symbols[i] = trim(match[1]);
        std::istringstream row(line);
        std::string token, last;
        while (row >> token) last = token;
        if (last.empty()) throw std::runtime_error("Bad species mass in dynG");
        type_masses[i] = number(last) / AMU_RY;
    }

    geometry.masses.resize(nat);
    geometry.r_bohr.resize(nat);
    geometry.symbols.resize(nat);
    for (int i = 0; i < nat; ++i) {
        if (!std::getline(input, line))
            throw std::runtime_error("Incomplete dynG atom block: " + path.string());
        std::istringstream row(line);
        int atom_index, type_index;
        Vec3 position{};
        row >> atom_index >> type_index >> position[0] >> position[1] >> position[2];
        if (!row || atom_index != i + 1 || type_index < 1 || type_index > ntyp ||
            !std::isfinite(position[0]) || !std::isfinite(position[1]) ||
            !std::isfinite(position[2]))
            throw std::runtime_error("Bad atom row in dynG: " + path.string());
        geometry.masses[i] = type_masses.at(type_index - 1);
        geometry.symbols[i] = type_symbols.at(type_index - 1);
        for (int k = 0; k < 3; ++k) geometry.r_bohr[i][k] = position[k] * alat;
    }
    return geometry;
}

} // namespace fdvib
