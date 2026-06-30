// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

#include "fdvib.hpp"

#include <iostream>
#include <stdexcept>

namespace fdvib {

namespace {

std::string execution_settings(const std::string &prefix) {
    return "displacement_angstrom = 0.01\n"
           "pw_command = pw.x\n"
           "prefix = " + prefix + "\n"
           "run_dynmat = true\n"
           "dynmat_command = dynmat.x\n";
}

std::string input_template(const std::string &type) {
    const std::string common =
        "# Change scf_input if your QE input uses another filename.\n"
        "scf_input = scf.in\n"
        "outdir = fdvib\n";

    if (type == "local") {
        return "# Replace selected_atoms with one-based QE atom indices.\n" + common +
               "system_type = local\n"
               "selected_atoms = 1\n" + execution_settings("system");
    }
    if (type == "gas") {
        return "# Set multiplicity to the molecular spin multiplicity.\n" + common +
               "system_type = gas\n"
               "selected_atoms = all\n"
               "multiplicity = 1\n" + execution_settings("molecule");
    }
    throw std::runtime_error("init type must be local or gas");
}

} // namespace

void initialize_input(const std::string &type, const fs::path &directory) {
    const auto contents = input_template(lower(type));
    const auto destination = directory / "fdvib.in";
    if (fs::exists(fs::symlink_status(destination)))
        throw std::runtime_error("Refuse to overwrite existing " + display_path(destination));

    write_text(destination, contents);
    std::cout << "Wrote " << display_path(destination) << '\n' << std::flush;
}

} // namespace fdvib
