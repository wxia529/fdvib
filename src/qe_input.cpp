// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

#include "fdvib.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace fdvib {

namespace {

std::string card_unit(const std::string &line, const std::string &card) {
    const std::regex re("^\\s*" + card + R"(\b\s*(?:[({]\s*([A-Za-z_]+)\s*[)}]|([A-Za-z_]+))?)",
                        std::regex::icase);
    std::smatch m;
    if (!std::regex_search(line, m, re)) return {};
    std::string unit = m[1].matched ? m[1].str() : (m[2].matched ? m[2].str() : "");
    return lower(unit);
}

double scalar_parameter(const std::string &text, const std::string &name) {
    std::smatch m;
    if (std::regex_search(text, m,
                          std::regex("\\b" + name + R"(\s*=\s*([-+0-9.EeDd]+))",
                                     std::regex::icase)))
        return number(m[1]);
    return 0.0;
}

double alat_angstrom(const std::string &text) {
    const double a = scalar_parameter(text, "A");
    if (a > 0.0) return a;
    const double celldm1 = scalar_parameter(text, R"(celldm\s*\(\s*1\s*\))");
    if (celldm1 > 0.0) return celldm1 * BOHR_TO_ANG;
    return 0.0;
}

Vec3 scaled(Vec3 v, double factor) {
    for (double &x : v) x *= factor;
    return v;
}

Vec3 crystal_to_cart(const Vec3 &frac, const std::array<Vec3, 3> &cell) {
    Vec3 out{0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            out[k] += frac[i] * cell[i][k];
    return out;
}

Vec3 cart_to_crystal_delta(const std::array<Vec3, 3> &cell, const Vec3 &cart) {
    const double a00 = cell[0][0], a01 = cell[1][0], a02 = cell[2][0];
    const double a10 = cell[0][1], a11 = cell[1][1], a12 = cell[2][1];
    const double a20 = cell[0][2], a21 = cell[1][2], a22 = cell[2][2];
    const double det = a00*(a11*a22 - a12*a21) - a01*(a10*a22 - a12*a20) +
                       a02*(a10*a21 - a11*a20);
    if (std::abs(det) < 1.0e-14) throw std::runtime_error("CELL_PARAMETERS matrix is singular");
    Vec3 out{};
    out[0] = ((a11*a22 - a12*a21)*cart[0] + (a02*a21 - a01*a22)*cart[1] +
              (a01*a12 - a02*a11)*cart[2]) / det;
    out[1] = ((a12*a20 - a10*a22)*cart[0] + (a00*a22 - a02*a20)*cart[1] +
              (a02*a10 - a00*a12)*cart[2]) / det;
    out[2] = ((a10*a21 - a11*a20)*cart[0] + (a01*a20 - a00*a21)*cart[1] +
              (a00*a11 - a01*a10)*cart[2]) / det;
    return out;
}

bool supported_position_unit(const std::string &unit) {
    return unit == "angstrom" || unit == "bohr" || unit == "alat" || unit == "crystal";
}

bool supported_cell_unit(const std::string &unit) {
    return unit == "angstrom" || unit == "bohr" || unit == "alat";
}

} // namespace

std::smatch search_or_throw(const std::string &text, const std::regex &re,
                            const std::string &what) {
    std::smatch m;
    if (!std::regex_search(text, m, re)) throw std::runtime_error("Cannot find " + what);
    return m;
}

QEInput parse_qe_input(const fs::path &p) {
    QEInput q;
    q.text = read_text(p);
    std::istringstream in(q.text);
    std::string line;
    while (std::getline(in, line)) {
        q.lines.push_back(line + "\n");
        q.clean_text += strip_comment(line) + "\n";
    }

    auto mnat = search_or_throw(q.clean_text, std::regex(R"(\bnat\s*=\s*(\d+))", std::regex::icase), "nat");
    auto mnt = search_or_throw(q.clean_text, std::regex(R"(\bntyp\s*=\s*(\d+))", std::regex::icase), "ntyp");
    q.nat = std::stoi(mnat[1]); q.ntyp = std::stoi(mnt[1]);
    if (q.nat <= 0 || q.ntyp <= 0) throw std::runtime_error("nat and ntyp must be positive");
    auto mib = search_or_throw(q.clean_text, std::regex(R"(\bibrav\s*=\s*([-+]?\d+))", std::regex::icase), "ibrav");
    if (std::stoi(mib[1]) != 0) throw std::runtime_error("FDVIB requires ibrav=0");
    if (!std::regex_search(q.clean_text, std::regex(R"(\bcalculation\s*=\s*['\"]scf['\"])", std::regex::icase)))
        throw std::runtime_error("scf.in must contain calculation='scf'");
    if (!std::regex_search(q.clean_text, std::regex(R"(\btprnfor\s*=\s*\.true\.)", std::regex::icase)))
        throw std::runtime_error("scf.in must contain tprnfor=.true.");
    const std::regex startingpot_re(R"(\bstartingpot\s*=\s*['\"]([^'\"]+)['\"])", std::regex::icase);
    for (std::sregex_iterator it(q.clean_text.begin(), q.clean_text.end(), startingpot_re), end;
         it != end; ++it)
        if (lower((*it)[1]) == "file")
            throw std::runtime_error("scf.in must not set startingpot='file'; FDVIB manages the reference density");
    std::smatch prefix;
    if (std::regex_search(q.clean_text, prefix,
                          std::regex(R"(\bprefix\s*=\s*['\"]([^'\"]+)['\"])", std::regex::icase)))
        q.prefix = prefix[1];
    if (q.prefix.empty() || fs::path(q.prefix).filename() != fs::path(q.prefix) ||
        q.prefix == "." || q.prefix == "..")
        throw std::runtime_error("QE prefix must be a non-empty filename prefix without directories");

    q.alat_angstrom = alat_angstrom(q.clean_text);
    int sp = -1;
    for (int i = 0; i < static_cast<int>(q.lines.size()); ++i) {
        const auto s = trim(q.lines[i]);
        if (std::regex_search(s, std::regex(R"(^ATOMIC_SPECIES\b)", std::regex::icase))) sp = i;
        if (std::regex_search(s, std::regex(R"(^ATOMIC_POSITIONS\b)", std::regex::icase))) {
            q.positions_unit = card_unit(s, "ATOMIC_POSITIONS");
            if (q.positions_unit.empty())
                throw std::runtime_error("ATOMIC_POSITIONS must specify units: angstrom, bohr, alat, or crystal");
            if (q.positions_unit == "crystal_sg")
                throw std::runtime_error("ATOMIC_POSITIONS crystal_sg is not supported; use explicit crystal coordinates");
            if (!supported_position_unit(q.positions_unit))
                throw std::runtime_error("Unsupported ATOMIC_POSITIONS unit: " + q.positions_unit);
            q.pos_header = i;
            q.pos_start = i + 1;
        }
        if (std::regex_search(s, std::regex(R"(^CELL_PARAMETERS\b)", std::regex::icase))) {
            q.cell_unit = card_unit(s, "CELL_PARAMETERS");
            if (q.cell_unit.empty())
                throw std::runtime_error("CELL_PARAMETERS must specify units: angstrom, bohr, or alat");
            if (!supported_cell_unit(q.cell_unit))
                throw std::runtime_error("Unsupported CELL_PARAMETERS unit: " + q.cell_unit);
            q.cell_header = i;
        }
    }
    if (sp < 0 || q.pos_start < 0 || q.cell_header < 0)
        throw std::runtime_error("Require ATOMIC_SPECIES, ATOMIC_POSITIONS, and CELL_PARAMETERS");
    const auto line_count = static_cast<int>(q.lines.size());
    if (sp + q.ntyp >= line_count)
        throw std::runtime_error("ATOMIC_SPECIES block is shorter than ntyp");
    if (q.pos_start + q.nat > line_count)
        throw std::runtime_error("ATOMIC_POSITIONS block is shorter than nat");
    if (q.cell_header + 3 >= line_count)
        throw std::runtime_error("CELL_PARAMETERS block is incomplete");
    if ((q.positions_unit == "alat" || q.cell_unit == "alat") && q.alat_angstrom <= 0.0)
        throw std::runtime_error("alat units require A or celldm(1) in &SYSTEM");

    for (int i = 0; i < q.ntyp; ++i) {
        std::istringstream ls(q.lines[sp + 1 + i]);
        Species x;
        if (!(ls >> x.symbol >> x.mass >> x.pseudo)) throw std::runtime_error("Bad ATOMIC_SPECIES line");
        if (!(x.mass > 0.0) || !std::isfinite(x.mass)) throw std::runtime_error("Atomic mass must be positive");
        q.species.push_back(x);
    }
    for (int i = 0; i < q.nat; ++i) {
        std::istringstream ls(q.lines[q.pos_start + i]);
        Atom a;
        if (!(ls >> a.symbol >> a.input_r[0] >> a.input_r[1] >> a.input_r[2])) throw std::runtime_error("Bad ATOMIC_POSITIONS line");
        if (!std::isfinite(a.input_r[0]) || !std::isfinite(a.input_r[1]) || !std::isfinite(a.input_r[2]))
            throw std::runtime_error("Atomic position must be finite");
        std::string x; while (ls >> x) a.extra.push_back(x);
        auto it = std::find_if(q.species.begin(), q.species.end(), [&](const Species &s){ return s.symbol == a.symbol; });
        if (it == q.species.end()) throw std::runtime_error("Unknown species " + a.symbol);
        a.type = static_cast<int>(it - q.species.begin()) + 1;
        q.atoms.push_back(a);
    }
    for (int i = 0; i < 3; ++i) {
        std::istringstream ls(q.lines[q.cell_header + 1 + i]);
        if (!(ls >> q.cell[i][0] >> q.cell[i][1] >> q.cell[i][2])) throw std::runtime_error("Bad CELL_PARAMETERS");
        if (!std::isfinite(q.cell[i][0]) || !std::isfinite(q.cell[i][1]) || !std::isfinite(q.cell[i][2]))
            throw std::runtime_error("CELL_PARAMETERS values must be finite");
    }
    const double cell_factor = q.cell_unit == "angstrom" ? 1.0 :
                               q.cell_unit == "bohr" ? BOHR_TO_ANG : q.alat_angstrom;
    for (auto &v : q.cell) v = scaled(v, cell_factor);
    for (auto &a : q.atoms) {
        if (q.positions_unit == "angstrom") a.r = a.input_r;
        else if (q.positions_unit == "bohr") a.r = scaled(a.input_r, BOHR_TO_ANG);
        else if (q.positions_unit == "alat") a.r = scaled(a.input_r, q.alat_angstrom);
        else if (q.positions_unit == "crystal") a.r = crystal_to_cart(a.input_r, q.cell);
    }
    return q;
}

std::string format_position(const Atom &a, const Vec3 &coord) {
    std::ostringstream o;
    o << std::left << std::setw(4) << a.symbol << std::right << std::fixed << std::setprecision(10)
      << ' ' << std::setw(18) << coord[0] << ' ' << std::setw(18) << coord[1]
      << ' ' << std::setw(18) << coord[2];
    for (const auto &x : a.extra) o << "  " << x;
    o << '\n';
    return o.str();
}

std::vector<std::string> attempt_lines(const QEInput &q, const std::string &outdir,
                                       const fs::path &source_dir, const fs::path &run_dir) {
    auto lines = q.lines;
    bool found = false;
    const std::regex outdir_re(R"(^\s*outdir\s*=)", std::regex::icase);
    const std::regex wfcdir_re(R"(^\s*wfcdir\s*=)", std::regex::icase);
    const std::regex pseudo_re(R"(^\s*pseudo_dir\s*=\s*(['\"])([^'\"]+)['\"])", std::regex::icase);
    for (auto &line : lines) {
        const auto indent = line.substr(0, line.find_first_not_of(" \t"));
        if (std::regex_search(line, outdir_re)) {
            line = indent + "outdir = '" + outdir + "',\n";
            found = true;
        } else if (std::regex_search(line, wfcdir_re)) {
            line = indent + "wfcdir = '" + outdir + "',\n";
        } else {
            std::smatch match;
            if (std::regex_search(line, match, pseudo_re)) {
                const auto original = match[2].str();
                fs::path value(original);
                const bool home_relative = original == "~" || original.rfind("~/", 0) == 0 ||
                                           original == "$HOME" || original.rfind("$HOME/", 0) == 0 ||
                                           original == "${HOME}" || original.rfind("${HOME}/", 0) == 0;
                if (home_relative) {
                    const char *home = std::getenv("HOME");
                    if (!home || !*home)
                        throw std::runtime_error("pseudo_dir uses the home directory but HOME is not set");
                    std::string suffix;
                    if (original.rfind("~/", 0) == 0) suffix = original.substr(2);
                    else if (original.rfind("$HOME/", 0) == 0) suffix = original.substr(6);
                    else if (original.rfind("${HOME}/", 0) == 0) suffix = original.substr(8);
                    value = (fs::path(home) / suffix).lexically_normal();
                } else if (!value.is_absolute() && !original.empty() &&
                           original.front() != '~' && original.front() != '$') {
                    const auto target = (source_dir / value).lexically_normal();
                    value = target.lexically_relative(run_dir.lexically_normal());
                    if (value.empty()) value = ".";
                }
                line = indent + "pseudo_dir = '" + value.generic_string() + "',\n";
            }
        }
    }
    if (!found) throw std::runtime_error("Cannot find outdir in scf.in");
    return lines;
}

std::string displaced_input(const QEInput &q, int atom, int axis, double shift,
                            const std::string &outdir, const fs::path &source_dir,
                            const fs::path &run_dir) {
    auto lines = attempt_lines(q, outdir, source_dir, run_dir);
    auto a = q.atoms.at(atom);
    Vec3 coord = a.input_r;
    if (q.positions_unit == "angstrom") {
        coord[axis] += shift;
    } else if (q.positions_unit == "bohr") {
        coord[axis] += shift / BOHR_TO_ANG;
    } else if (q.positions_unit == "alat") {
        coord[axis] += shift / q.alat_angstrom;
    } else if (q.positions_unit == "crystal") {
        Vec3 cart{0.0, 0.0, 0.0};
        cart[axis] = shift;
        const auto delta = cart_to_crystal_delta(q.cell, cart);
        for (int k = 0; k < 3; ++k) coord[k] += delta[k];
    }
    lines[q.pos_start + atom] = format_position(a, coord);
    bool disk_io_found = false;
    const std::regex disk_io_re(R"(^\s*disk_io\s*=)", std::regex::icase);
    for (auto &line : lines) if (std::regex_search(line, disk_io_re)) {
        const auto indent = line.substr(0, line.find_first_not_of(" \t"));
        line = indent + "disk_io = 'nowf',\n";
        disk_io_found = true;
    }
    if (!disk_io_found) {
        int control = -1;
        for (int i = 0; i < static_cast<int>(lines.size()); ++i)
            if (std::regex_search(trim(lines[i]), std::regex(R"(^&CONTROL\b)", std::regex::icase))) {
                control = i;
                break;
            }
        if (control < 0) throw std::runtime_error("Cannot find &CONTROL in scf.in");
        int end = control + 1;
        while (end < static_cast<int>(lines.size()) && trim(strip_comment(lines[end])) != "/") ++end;
        if (end == static_cast<int>(lines.size())) throw std::runtime_error("Unterminated &CONTROL namelist");
        lines.insert(lines.begin() + end, "  disk_io = 'nowf',\n");
    }
    bool startingpot_found = false;
    std::regex spre(R"(^\s*startingpot\s*=)", std::regex::icase);
    for (auto &l : lines) if (std::regex_search(l, spre)) {
        const auto indent = l.substr(0, l.find_first_not_of(" \t"));
        l = indent + "startingpot = 'file',\n";
        startingpot_found = true;
    }
    if (!startingpot_found) {
        int electrons = -1;
        for (int i = 0; i < static_cast<int>(lines.size()); ++i)
            if (std::regex_search(trim(lines[i]), std::regex(R"(^&ELECTRONS\b)", std::regex::icase))) {
                electrons = i;
                break;
            }
        if (electrons >= 0) {
            int end = electrons + 1;
            while (end < static_cast<int>(lines.size()) && trim(strip_comment(lines[end])) != "/") ++end;
            if (end == static_cast<int>(lines.size())) throw std::runtime_error("Unterminated &ELECTRONS namelist");
            lines.insert(lines.begin() + end, "  startingpot = 'file',\n");
        } else {
            int first_card = std::min({q.pos_header, q.cell_header});
            for (int i = 0; i < static_cast<int>(lines.size()); ++i)
                if (std::regex_search(trim(lines[i]), std::regex(R"(^ATOMIC_SPECIES\b)", std::regex::icase)))
                    first_card = std::min(first_card, i);
            lines.insert(lines.begin() + first_card, "&ELECTRONS\n  startingpot = 'file',\n/\n");
        }
    }
    return std::accumulate(lines.begin(), lines.end(), std::string{});
}

std::string reference_input(const QEInput &q, const std::string &outdir,
                            const fs::path &source_dir, const fs::path &run_dir) {
    const auto lines = attempt_lines(q, outdir, source_dir, run_dir);
    return std::accumulate(lines.begin(), lines.end(), std::string{});
}

} // namespace fdvib
