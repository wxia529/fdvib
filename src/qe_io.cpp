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

    int sp = -1;
    for (int i = 0; i < static_cast<int>(q.lines.size()); ++i) {
        const auto s = trim(q.lines[i]);
        if (std::regex_search(s, std::regex(R"(^ATOMIC_SPECIES\b)", std::regex::icase))) sp = i;
        if (std::regex_search(s, std::regex(R"(^ATOMIC_POSITIONS\s*[({]?\s*angstrom\s*[)}]?)", std::regex::icase))) {
            q.pos_header = i; q.pos_start = i + 1;
        }
        if (std::regex_search(s, std::regex(R"(^CELL_PARAMETERS\s*[({]?\s*angstrom\s*[)}]?)", std::regex::icase)))
            q.cell_header = i;
    }
    if (sp < 0 || q.pos_start <= 0 || q.cell_header <= 0)
        throw std::runtime_error("Require ATOMIC_SPECIES, ATOMIC_POSITIONS angstrom, CELL_PARAMETERS angstrom");

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
        if (!(ls >> a.symbol >> a.r[0] >> a.r[1] >> a.r[2])) throw std::runtime_error("Bad ATOMIC_POSITIONS line");
        std::string x; while (ls >> x) a.extra.push_back(x);
        auto it = std::find_if(q.species.begin(), q.species.end(), [&](const Species &s){ return s.symbol == a.symbol; });
        if (it == q.species.end()) throw std::runtime_error("Unknown species " + a.symbol);
        a.type = static_cast<int>(it - q.species.begin()) + 1;
        q.atoms.push_back(a);
    }
    for (int i = 0; i < 3; ++i) {
        std::istringstream ls(q.lines[q.cell_header + 1 + i]);
        if (!(ls >> q.cell[i][0] >> q.cell[i][1] >> q.cell[i][2])) throw std::runtime_error("Bad CELL_PARAMETERS");
    }
    return q;
}

std::string format_position(const Atom &a) {
    std::ostringstream o;
    o << std::left << std::setw(4) << a.symbol << std::right << std::fixed << std::setprecision(10)
      << ' ' << std::setw(18) << a.r[0] << ' ' << std::setw(18) << a.r[1]
      << ' ' << std::setw(18) << a.r[2];
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
    a.r[axis] += shift;
    lines[q.pos_start + atom] = format_position(a);
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
