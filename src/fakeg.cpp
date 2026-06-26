#include "fdvib.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace fdvib {

namespace {

int atomic_number(const std::string &label) {
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
    const auto find_symbol = [&](const std::string &symbol) -> int {
        for (std::size_t z = 1; z < symbols.size(); ++z)
            if (symbol == symbols[z]) return static_cast<int>(z);
        return 0;
    };
    if (const int z = find_symbol(label); z) return z;
    if (label.size() >= 2)
        if (const int z = find_symbol(label.substr(0, 2)); z) return z;
    if (!label.empty())
        if (const int z = find_symbol(label.substr(0, 1)); z) return z;
    throw std::runtime_error("Cannot map species label to atomic number for fake Gaussian export: " + label);
}

void write_standard_orientation(std::ostringstream &out, const DynGeometry &geometry) {
    out << "                         Standard orientation:\n"
        << " ---------------------------------------------------------------------\n"
        << " Center     Atomic      Atomic             Coordinates (Angstroms)\n"
        << " Number     Number       Type             X           Y           Z\n"
        << " ---------------------------------------------------------------------\n";
    for (std::size_t i = 0; i < geometry.symbols.size(); ++i) {
        out << std::setw(7) << i + 1
            << std::setw(11) << atomic_number(geometry.symbols[i])
            << std::setw(12) << 0
            << std::fixed << std::setprecision(6)
            << std::setw(16) << geometry.r_bohr[i][0] * BOHR_TO_ANG
            << std::setw(12) << geometry.r_bohr[i][1] * BOHR_TO_ANG
            << std::setw(12) << geometry.r_bohr[i][2] * BOHR_TO_ANG << '\n';
    }
    out << " ---------------------------------------------------------------------\n\n";
}

} // namespace

void fakeg(const fs::path &results) {
    const auto files = result_files(results, "Fake Gaussian export");
    const auto geometry = read_dyn_geometry(files.dyn);
    const auto modes = parse_modes(files.freq, static_cast<int>(geometry.masses.size()));
    const auto destination = results / (files.dyn.stem().string() + "_fake.out");
    if (fs::exists(destination)) throw std::runtime_error("Refuse to overwrite " + destination.string());

    std::ostringstream out;
    out << " 0 basis functions\n"
        << " 0 alpha electrons\n"
        << " 0 beta electrons\n"
        << "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n"
        << "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n\n";
    write_standard_orientation(out, geometry);
    out << " SCF Done:         0.000000000\n\n"
        << " Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n"
        << " activities (A**4/AMU), depolarization ratios for plane and unpolarized\n"
        << " incident light, reduced masses (AMU), force constants (mDyne/A),\n"
        << " and normal coordinates:\n";

    for (std::size_t first = 0; first < modes.size(); first += 3) {
        const std::size_t count = std::min<std::size_t>(3, modes.size() - first);
        out << "                 ";
        for (std::size_t j = 0; j < count; ++j) out << std::setw(23) << first + j + 1;
        out << "\n                 ";
        for (std::size_t j = 0; j < count; ++j) out << std::setw(23) << "A";
        out << "\n Frequencies --";
        for (std::size_t j = 0; j < count; ++j)
            out << std::fixed << std::setprecision(4) << std::setw(12) << modes[first + j].freq << "           ";
        out << "\n IR Inten    --";
        for (std::size_t j = 0; j < count; ++j)
            out << std::fixed << std::setprecision(4) << std::setw(12) << 0.0 << "           ";
        out << "\n  Atom  AN      X      Y      Z  ";
        for (std::size_t j = 1; j < count; ++j) out << "      X      Y      Z  ";
        out << '\n';
        for (std::size_t atom = 0; atom < geometry.symbols.size(); ++atom) {
            out << std::setw(6) << atom + 1 << std::setw(4) << atomic_number(geometry.symbols[atom]);
            for (std::size_t j = 0; j < count; ++j) {
                const auto &d = modes[first + j].displacement.at(atom);
                out << std::fixed << std::setprecision(2)
                    << std::setw(9) << d[0] << std::setw(7) << d[1] << std::setw(7) << d[2];
            }
            out << '\n';
        }
    }
    out << "\n Temperature     0.000 Kelvin.  Pressure   0.00000 Atm.\n"
        << " Zero-point correction=                           0.000000 Hartree\n"
        << " Thermal correction to Energy=                    0.000000\n"
        << " Thermal correction to Enthalpy=                  0.000000\n"
        << " Thermal correction to Gibbs Free Energy=         0.000000\n"
        << " Normal termination of Gaussian\n";

    write_text(destination, out.str());
    std::cout << "Wrote " << destination << "\n";
}

} // namespace fdvib
