#include "fdvib.hpp"

#include <algorithm>
#include <iomanip>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace fdvib {

void validate_qe_output(const fs::path &output) {
    const auto text = read_text(output);
    if (text.find("JOB DONE") == std::string::npos)
        throw std::runtime_error("JOB DONE not found: " + output.string());
    if (std::regex_search(text, std::regex(R"(convergence\s+NOT\s+achieved)", std::regex::icase)))
        throw std::runtime_error("SCF convergence was not achieved: " + output.string());
    if (std::regex_search(text, std::regex(R"(Error\s+in\s+routine)", std::regex::icase)))
        throw std::runtime_error("QE reported Error in routine: " + output.string());
}

std::vector<Vec3> parse_forces(const fs::path &output, int nat) {
    validate_qe_output(output);
    const auto text = read_text(output);
    const auto pos = text.rfind("Forces acting on atoms");
    if (pos == std::string::npos) throw std::runtime_error("Force block not found: " + output.string());
    std::istringstream in(text.substr(pos));
    std::regex re(R"(^\s*atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+))", std::regex::icase);
    std::vector<Vec3> f(nat);
    std::vector<bool> got(nat, false);
    std::string line;
    std::smatch m;
    bool started = false;
    int count = 0;
    while (std::getline(in, line)) {
        if (std::regex_search(line, m, re)) {
            started = true;
            const int i = std::stoi(m[1]) - 1;
            if (i < 0 || i >= nat || got[i])
                throw std::runtime_error("Invalid/duplicate atom in force block: " + output.string());
            f[i] = {number(m[2]), number(m[3]), number(m[4])};
            got[i] = true;
            if (++count == nat) return f;
        } else if (started) {
            break;
        }
    }
    throw std::runtime_error("Incomplete force block: " + output.string());
}

void write_forces(const fs::path &p, const std::vector<Vec3> &f, const fs::path &source) {
    std::ostringstream o;
    o << "# source: " << source.filename().string()
      << "\n# nat: " << f.size()
      << "\n# units: Ry/Bohr\n# atom fx fy fz\n";
    o << std::scientific << std::setprecision(15);
    for (std::size_t i = 0; i < f.size(); ++i)
        o << i + 1 << ' ' << f[i][0] << ' ' << f[i][1] << ' ' << f[i][2] << '\n';
    write_text(p, o.str());
}

std::vector<Vec3> read_forces(const fs::path &p, int nat) {
    std::istringstream in(read_text(p));
    std::vector<Vec3> f(nat);
    std::vector<bool> got(nat, false);
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line.front() == '#') continue;
        std::istringstream ls(line);
        int i;
        Vec3 v;
        if (!(ls >> i >> v[0] >> v[1] >> v[2]) || i < 1 || i > nat)
            throw std::runtime_error("Bad forces.dat: " + p.string());
        f[i - 1] = v;
        got[i - 1] = true;
    }
    if (std::count(got.begin(), got.end(), true) != nat)
        throw std::runtime_error("Incomplete forces.dat: " + p.string());
    return f;
}

double read_total_energy_hartree(const fs::path &output) {
    const auto text = read_text(output);
    const std::regex pattern(R"(!\s+total\s+energy\s*=\s*([-+0-9.EeDd]+)\s+Ry\b)",
                             std::regex::icase);
    bool found = false;
    double energy_ry = 0.0;
    for (std::sregex_iterator it(text.begin(), text.end(), pattern), end; it != end; ++it) {
        energy_ry = number((*it)[1]);
        found = true;
    }
    if (!found) throw std::runtime_error("Cannot find converged QE total energy in " + output.string());
    return energy_ry / 2.0;
}

} // namespace fdvib
