#include "fdvib.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace fdvib {

namespace {

void require_no_extra(std::istringstream &in) {
    std::string extra;
    if (in >> extra) throw std::runtime_error("Extra field in generated SHM data line");
}

void validate_shm(const fs::path &path) {
    std::istringstream input(read_text(path));
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) throw std::runtime_error("Generated SHM contains an empty line");
        lines.push_back(line);
    }
    if (lines.empty() || read_text(path).back() != '\n') throw std::runtime_error("Generated SHM must end with a newline");
    std::vector<std::size_t> tags;
    const std::vector<std::string> names{"*E", "*wavenum", "*atoms", "*elevel"};
    for (const auto &name : names) {
        std::vector<std::size_t> found;
        for (std::size_t i = 0; i < lines.size(); ++i) if (lines[i].rfind(name, 0) == 0) found.push_back(i);
        if (found.size() != 1) throw std::runtime_error("Generated SHM must contain exactly one " + name + " tag");
        tags.push_back(found.front());
    }
    if (!(tags[0] + 1 < tags[1] && tags[1] < tags[2] && tags[2] + 1 < tags[3]))
        throw std::runtime_error("Generated SHM sections are missing data or out of order");
    if (tags[1] != tags[0] + 2) throw std::runtime_error("Generated SHM *E section must contain one value");
    {
        std::istringstream row(lines[tags[0]+1]); double value;
        if (!(row >> value) || !std::isfinite(value)) throw std::runtime_error("Invalid SHM electronic energy");
        require_no_extra(row);
    }
    for (std::size_t i = tags[1]+1; i < tags[2]; ++i) {
        std::istringstream row(lines[i]); double value;
        if (!(row >> value) || !std::isfinite(value)) throw std::runtime_error("Invalid SHM wavenumber");
        require_no_extra(row);
    }
    for (std::size_t i = tags[2]+1; i < tags[3]; ++i) {
        std::istringstream row(lines[i]); std::string symbol; double mass,x,y,z;
        if (!(row >> symbol >> mass >> x >> y >> z) || !is_standard_element(symbol) ||
            !(mass > 0.0) || !std::isfinite(mass) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            throw std::runtime_error("Invalid SHM atom row");
        require_no_extra(row);
    }
    bool has_ground = false;
    for (std::size_t i = tags[3]+1; i < lines.size(); ++i) {
        std::istringstream row(lines[i]); double energy; int degeneracy;
        if (!(row >> energy >> degeneracy) || !std::isfinite(energy) || energy < 0.0 || degeneracy <= 0)
            throw std::runtime_error("Invalid SHM electronic-level row");
        require_no_extra(row);
        if (energy == 0.0) has_ground = true;
    }
    if (tags[3]+1 == lines.size() || !has_ground) throw std::runtime_error("Generated SHM requires a ground electronic level");
}

} // namespace

void shm(const fs::path &results) {
    const auto files = result_files(results, "SHM export");
    const auto geometry = read_dyn_geometry(files.dyn);
    const auto modes = parse_modes(files.freq, static_cast<int>(geometry.masses.size()));
    const auto metadata = result_metadata(results, false);
    const double energy = metadata.electronic_energy_hartree;
    const int multiplicity = metadata.multiplicity;
    const auto selected = select_shm_modes(metadata, geometry, modes);

    double smallest_retained = 0.0;
    if (!selected.modes.empty()) {
        smallest_retained = std::numeric_limits<double>::infinity();
        for (const auto *mode : selected.modes)
            smallest_retained = std::min(smallest_retained, std::abs(mode->freq));
    }
    std::vector<std::string> output_symbols;
    output_symbols.reserve(geometry.symbols.size());
    for (std::size_t i = 0; i < geometry.symbols.size(); ++i) {
        output_symbols.push_back(standard_element_symbol(geometry.symbols[i], "SHM export"));
        if (!(geometry.masses[i] > 0.0)) throw std::runtime_error("SHM export requires positive atomic masses");
    }

    std::ostringstream output;
    output << "*E  //Electronic energy (a.u.)\n" << std::uppercase << std::scientific << std::setprecision(15) << energy << '\n'
           << "*wavenum  //Wavenumbers (cm-1). Negative value means imaginary frequency\n";
    for (const auto *mode : selected.modes) output << std::fixed << std::setprecision(10) << mode->freq << '\n';
    output << "*atoms  //Information of all atoms: Name, mass (amu), X, Y, Z (Angstrom)\n";
    for (std::size_t i = 0; i < geometry.symbols.size(); ++i)
        output << output_symbols[i] << ' ' << std::fixed << std::setprecision(10) << geometry.masses[i] << ' '
               << geometry.r_bohr[i][0]*BOHR_TO_ANG << ' ' << geometry.r_bohr[i][1]*BOHR_TO_ANG << ' '
               << geometry.r_bohr[i][2]*BOHR_TO_ANG << '\n';
    output << "*elevel  //Energy (eV) and degeneracy of electronic energy levels\n"
           << "0.000000 " << multiplicity << '\n';

    const auto destination = results / (files.dyn.stem().string() + ".shm");
    const auto temporary = fs::path(destination.string() + ".tmp");
    if (fs::exists(temporary)) throw std::runtime_error("Stale SHM temporary file exists: " + temporary.string());
    write_text(temporary, output.str());
    validate_shm(temporary);
    fs::rename(temporary, destination);
    std::cout << "Wrote " << destination << "\n"
              << "SHM mode selection: " << selected.classification << ", retained " << selected.modes.size()
              << ", removed " << modes.size()-selected.modes.size()
              << ", largest removed |frequency|=" << selected.largest_removed
              << ", smallest retained |frequency|=" << smallest_retained << " cm^-1\n";
    if (metadata.mode_selection == "local")
        std::cout << "Run Shermo with: Shermo " << destination << " -imode 1 -PGlabel C1\n";
    else
        std::cout << "Run Shermo with: Shermo " << destination << "\n";
}

} // namespace fdvib
