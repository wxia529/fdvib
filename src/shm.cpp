#include "fdvib.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace fdvib {

namespace {

const std::set<std::string> &element_symbols() {
    static const std::set<std::string> symbols = {
        "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
        "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
        "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
        "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
        "Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    };
    return symbols;
}

std::string canonical_element(const std::string &label) {
    if (element_symbols().count(label)) return label;
    if (label.size() >= 2) {
        const auto two = label.substr(0, 2);
        if (element_symbols().count(two)) return two;
    }
    if (!label.empty()) {
        const auto one = label.substr(0, 1);
        if (element_symbols().count(one)) return one;
    }
    throw std::runtime_error("Cannot map QE species label to a standard element for SHM export: " + label);
}

std::array<double, 3> inertia_eigenvalues(const DynGeometry &g) {
    const double total_mass = std::accumulate(g.masses.begin(), g.masses.end(), 0.0);
    if (!(total_mass > 0.0)) throw std::runtime_error("Invalid total mass for SHM export");
    Vec3 center{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < g.masses.size(); ++i)
        for (int k = 0; k < 3; ++k) center[k] += g.masses[i] * g.r_bohr[i][k] / total_mass;
    Vec3 diagonal{0.0, 0.0, 0.0}, off{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < g.masses.size(); ++i) {
        const double x = g.r_bohr[i][0] - center[0];
        const double y = g.r_bohr[i][1] - center[1];
        const double z = g.r_bohr[i][2] - center[2];
        const double mass = g.masses[i];
        diagonal[0] += mass * (y*y + z*z);
        diagonal[1] += mass * (x*x + z*z);
        diagonal[2] += mass * (x*x + y*y);
        off[0] -= mass*x*y; off[1] -= mass*x*z; off[2] -= mass*y*z;
    }
    const double p1 = off[0]*off[0] + off[1]*off[1] + off[2]*off[2];
    if (p1 == 0.0) {
        std::array<double, 3> values{diagonal[0], diagonal[1], diagonal[2]};
        std::sort(values.begin(), values.end());
        return values;
    }
    const double q = (diagonal[0] + diagonal[1] + diagonal[2]) / 3.0;
    const double p2 = (diagonal[0]-q)*(diagonal[0]-q) + (diagonal[1]-q)*(diagonal[1]-q) +
                      (diagonal[2]-q)*(diagonal[2]-q) + 2.0*p1;
    const double p = std::sqrt(p2 / 6.0);
    const double a00=(diagonal[0]-q)/p, a11=(diagonal[1]-q)/p, a22=(diagonal[2]-q)/p;
    const double a01=off[0]/p, a02=off[1]/p, a12=off[2]/p;
    const double determinant = a00*a11*a22 + 2*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
    const double phi = std::acos(std::clamp(determinant/2.0, -1.0, 1.0)) / 3.0;
    const double e3=q+2*p*std::cos(phi), e1=q+2*p*std::cos(phi+2*PI/3), e2=3*q-e1-e3;
    std::array<double, 3> values{e1,e2,e3};
    std::sort(values.begin(), values.end());
    return values;
}

std::vector<double> retain_largest_modes(const std::vector<Mode> &modes, std::size_t keep,
                                         double &largest_removed) {
    if (keep > modes.size()) throw std::runtime_error("Requested SHM vibration count exceeds available modes");
    std::vector<std::size_t> order(modes.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return std::abs(modes[a].freq) < std::abs(modes[b].freq);
    });
    std::vector<bool> removed(modes.size(), false);
    largest_removed = 0.0;
    for (std::size_t i = 0; i < modes.size() - keep; ++i) {
        removed[order[i]] = true;
        largest_removed = std::max(largest_removed, std::abs(modes[order[i]].freq));
    }
    std::vector<double> frequencies;
    frequencies.reserve(keep);
    for (std::size_t i = 0; i < modes.size(); ++i)
        if (!removed[i]) frequencies.push_back(modes[i].freq);
    return frequencies;
}

std::vector<int> metadata_selected_atoms(const ResultMetadata &metadata, std::size_t nat) {
    std::vector<int> selected;
    if (metadata.selected_atoms == "all") {
        selected.resize(nat);
        std::iota(selected.begin(), selected.end(), 1);
    } else {
        selected = integer_list(metadata.selected_atoms);
    }
    if (selected.empty()) throw std::runtime_error("Local SHM export requires selected_atoms in metadata.dat");
    std::set<int> unique;
    for (const int atom : selected)
        if (atom < 1 || static_cast<std::size_t>(atom) > nat || !unique.insert(atom).second)
            throw std::runtime_error("Invalid/duplicate selected atom in metadata.dat");
    return selected;
}

struct ShmFrequencies {
    std::vector<double> values;
    std::string classification;
    double largest_removed{};
};

ShmFrequencies select_shm_frequencies(const ResultMetadata &metadata,
                                      const DynGeometry &geometry,
                                      const std::vector<Mode> &modes) {
    ShmFrequencies selected;
    const auto nat = geometry.masses.size();
    if (metadata.mode_selection == "all") {
        selected.classification = "all";
        selected.values.reserve(modes.size());
        for (const auto &mode : modes)
            if (mode.freq != 0.0) selected.values.push_back(mode.freq);
        return selected;
    }
    std::size_t keep = 0;
    if (metadata.mode_selection == "gas") {
        if (nat == 1) {
            keep = 0;
            selected.classification = "atom";
        } else {
            const auto inertia = inertia_eigenvalues(geometry);
            const bool linear = std::any_of(inertia.begin(), inertia.end(), [](double value) { return value < 0.001; });
            keep = 3*nat - (linear ? 5 : 6);
            selected.classification = linear ? "linear" : "nonlinear";
        }
    } else if (metadata.mode_selection == "local") {
        keep = 3 * metadata_selected_atoms(metadata, nat).size();
        selected.classification = "local";
    } else {
        throw std::runtime_error("mode_selection must be all, gas, or local in metadata.dat");
    }
    selected.values = retain_largest_modes(modes, keep, selected.largest_removed);
    return selected;
}

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
        if (!(row >> symbol >> mass >> x >> y >> z) || !element_symbols().count(symbol) ||
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
    const auto selected = select_shm_frequencies(metadata, geometry, modes);

    double smallest_retained = 0.0;
    if (!selected.values.empty()) {
        smallest_retained = std::numeric_limits<double>::infinity();
        for (const double frequency : selected.values)
            smallest_retained = std::min(smallest_retained, std::abs(frequency));
    }
    std::vector<std::string> output_symbols;
    output_symbols.reserve(geometry.symbols.size());
    for (std::size_t i = 0; i < geometry.symbols.size(); ++i) {
        output_symbols.push_back(canonical_element(geometry.symbols[i]));
        if (!(geometry.masses[i] > 0.0)) throw std::runtime_error("SHM export requires positive atomic masses");
    }

    std::ostringstream output;
    output << "*E  //Electronic energy (a.u.)\n" << std::uppercase << std::scientific << std::setprecision(15) << energy << '\n'
           << "*wavenum  //Wavenumbers (cm-1). Negative value means imaginary frequency\n";
    for (const double frequency : selected.values) output << std::fixed << std::setprecision(10) << frequency << '\n';
    output << "*atoms  //Information of all atoms: Name, mass (amu), X, Y, Z (Angstrom)\n";
    for (std::size_t i = 0; i < geometry.symbols.size(); ++i)
        output << output_symbols[i] << ' ' << std::fixed << std::setprecision(10) << geometry.masses[i] << ' '
               << geometry.r_bohr[i][0]*BOHR_TO_ANG << ' ' << geometry.r_bohr[i][1]*BOHR_TO_ANG << ' '
               << geometry.r_bohr[i][2]*BOHR_TO_ANG << '\n';
    output << "*elevel  //Energy (eV) and degeneracy of electronic energy levels\n"
           << "0.000000 " << multiplicity << '\n';

    const auto destination = results / (files.dyn.stem().string() + ".shm");
    const auto temporary = fs::path(destination.string() + ".tmp");
    if (fs::exists(destination)) throw std::runtime_error("Refuse to overwrite " + destination.string());
    if (fs::exists(temporary)) throw std::runtime_error("Stale SHM temporary file exists: " + temporary.string());
    write_text(temporary, output.str());
    validate_shm(temporary);
    fs::rename(temporary, destination);
    std::cout << "Wrote " << destination << "\n"
              << "SHM mode selection: " << selected.classification << ", retained " << selected.values.size()
              << ", removed " << modes.size()-selected.values.size()
              << ", largest removed |frequency|=" << selected.largest_removed
              << ", smallest retained |frequency|=" << smallest_retained << " cm^-1\n";
    if (metadata.mode_selection == "local")
        std::cout << "Run Shermo with: Shermo " << destination << " -imode 1 -PGlabel C1\n";
    else
        std::cout << "Run Shermo with: Shermo " << destination << "\n";
}

} // namespace fdvib
