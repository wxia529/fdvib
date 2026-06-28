// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

#include "fdvib.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

namespace fdvib {

namespace {

void write_molden(const fs::path &path, const DynGeometry &geometry,
                  const std::vector<Mode> &modes, const ModeSelection &selected,
                  bool include_cell) {
    std::vector<std::string> symbols;
    symbols.reserve(geometry.symbols.size());
    for (const auto &label : geometry.symbols)
        symbols.push_back(standard_element_symbol(label, "Molden export"));

    std::ostringstream output;
    output << "[Molden Format]\n" << std::fixed << std::setprecision(8);
    if (include_cell) {
        output << "[Cell]\n";
        for (const auto &vector : geometry.cell_bohr)
            output << std::setw(15) << vector[0] * BOHR_TO_ANG
                   << std::setw(15) << vector[1] * BOHR_TO_ANG
                   << std::setw(15) << vector[2] * BOHR_TO_ANG << '\n';
    }
    output << "[Atoms] AU\n";
    for (std::size_t i = 0; i < symbols.size(); ++i)
        output << std::setw(6) << symbols[i] << std::setw(8) << i + 1
               << std::setw(8) << atomic_number_from_label(symbols[i], "Molden export")
               << std::setw(18) << geometry.r_bohr[i][0]
               << std::setw(18) << geometry.r_bohr[i][1]
               << std::setw(18) << geometry.r_bohr[i][2] << '\n';
    output << "[FREQ]\n";
    for (const auto index : selected.indices) output << modes.at(index).freq << '\n';
    output << "[FR-COORD]\n";
    for (std::size_t i = 0; i < symbols.size(); ++i)
        output << std::setw(6) << symbols[i]
               << std::setw(18) << geometry.r_bohr[i][0]
               << std::setw(18) << geometry.r_bohr[i][1]
               << std::setw(18) << geometry.r_bohr[i][2] << '\n';
    output << "[FR-NORM-COORD]\n" << std::setprecision(10);
    for (std::size_t i = 0; i < selected.indices.size(); ++i) {
        output << " vibration " << i + 1 << '\n';
        for (const auto &displacement : modes.at(selected.indices[i]).displacement)
            output << std::setw(15) << displacement[0]
                   << std::setw(15) << displacement[1]
                   << std::setw(15) << displacement[2] << '\n';
    }
    write_text(path, output.str());
    std::cout << "Wrote " << selected.indices.size() << " Molden modes to "
              << display_path(path) << '\n';
}

} // namespace

void modes(const fs::path &results) {
    const auto files = result_files(results, "Normal-mode analysis");
    const auto geometry = read_qe_dyn_geometry(files.dyn);
    const auto parsed = read_qe_dynmat_modes(files.freq, static_cast<int>(geometry.masses.size()));
    const auto metadata = result_metadata(results, false);
    const auto selected = select_molden_modes(metadata, geometry, parsed);
    const auto output = results / (files.dyn.stem().string() + ".mol");
    write_molden(output, geometry, parsed, selected, metadata.mode_selection != "gas");
}

} // namespace fdvib
