#include "fdvib.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace fdvib {

std::string trim(const std::string &s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::string unquote(std::string s) {
    s = trim(s);
    if (s.size() >= 2 && ((s.front() == '\'' && s.back() == '\'') ||
                          (s.front() == '"' && s.back() == '"'))) {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

std::string strip_comment(const std::string &line) {
    bool single = false, dbl = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '\'' && !dbl) single = !single;
        else if (line[i] == '"' && !single) dbl = !dbl;
        else if (line[i] == '!' && !single && !dbl) return line.substr(0, i);
    }
    return line;
}

double number(std::string s) {
    s = unquote(trim(s));
    std::replace(s.begin(), s.end(), 'D', 'E');
    std::replace(s.begin(), s.end(), 'd', 'e');
    std::size_t used = 0;
    const double x = std::stod(s, &used);
    if (used != s.size()) throw std::runtime_error("Bad number: " + s);
    if (!std::isfinite(x)) throw std::runtime_error("Non-finite number: " + s);
    return x;
}

std::string read_text(const fs::path &p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot read " + p.string());
    std::ostringstream out;
    out << f.rdbuf();
    return out.str();
}

void write_text(const fs::path &p, const std::string &s) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write " + p.string());
    f << s;
    if (!f) throw std::runtime_error("Write failed: " + p.string());
}

Config Config::load(const fs::path &p) {
    Config c;
    std::istringstream in(read_text(p));
    std::string line;
    while (std::getline(in, line)) {
        line = trim(strip_comment(line));
        if (line.empty() || line.front() == '&' || line == "/") continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        auto key = lower(trim(line.substr(0, eq)));
        auto val = trim(line.substr(eq + 1));
        if (!val.empty() && val.back() == ',') val.pop_back();
        c.v[key] = trim(val);
    }
    return c;
}

bool Config::has(const std::string &k) const { return v.count(lower(k)) != 0; }
std::string Config::get(const std::string &k, const std::string &d) const {
    const auto it = v.find(lower(k));
    return it == v.end() ? d : unquote(it->second);
}
double Config::real(const std::string &k, double d) const {
    return has(k) ? number(v.at(lower(k))) : d;
}
int Config::integer(const std::string &k, int d) const {
    if (!has(k)) return d;
    const double x = number(v.at(lower(k)));
    const double rounded = std::round(x);
    if (std::abs(x - rounded) > 1.0e-10)
        throw std::runtime_error("Expected integer for " + k + ": " + v.at(lower(k)));
    return static_cast<int>(rounded);
}

void Config::require_only(const std::set<std::string> &allowed,
                          const std::string &context) const {
    for (const auto &[key, value] : v) {
        (void)value;
        if (!allowed.count(key))
            throw std::runtime_error("Unknown parameter in " + context + ": " + key);
    }
}

std::vector<int> integer_list(std::string s) {
    for (char &c : s) if (c == ',' || c == ';') c = ' ';
    std::istringstream in(s);
    std::vector<int> out;
    int x;
    while (in >> x) out.push_back(x);
    return out;
}

Settings settings(const fs::path &config_path, const fs::path &root_override) {
    Settings s;
    s.config_path = fs::absolute(config_path);
    s.root = root_override.empty() ? s.config_path.parent_path() : fs::absolute(root_override);
    const auto c = Config::load(s.config_path);
    c.require_only({"scf_input", "workdir", "system_type", "selected_atoms",
                    "displacement_angstrom", "multiplicity", "pw_command",
                    "mpi_command", "output_prefix"}, "fdvib.in");
    s.scf_input = s.root / c.get("scf_input", "scf.in");
    s.workdir = s.root / c.get("workdir", "fdvib");
    s.system_type = lower(c.get("system_type", "local"));
    s.displacement = c.real("displacement_angstrom", 0.01);
    s.multiplicity = c.integer("multiplicity", 1);
    s.multiplicity_explicit = c.has("multiplicity");
    s.pw_command = c.get("pw_command", "pw.x");
    s.mpi_command = c.get("mpi_command", "");
    s.output_prefix = c.get("output_prefix", "system");
    if (s.displacement <= 0) throw std::runtime_error("displacement_angstrom must be positive");
    if (s.system_type != "gas" && s.system_type != "local")
        throw std::runtime_error("system_type must be gas or local");
    const auto atoms = lower(c.get("selected_atoms", ""));
    s.selected_all = atoms == "all";
    if (s.selected_all) s.selected.clear();
    else s.selected = integer_list(atoms);
    return s;
}

std::string job_name(int atom1, int axis, int sign) {
    static const char *xyz[] = {"x", "y", "z"};
    std::ostringstream o;
    o << "disp_" << std::setfill('0') << std::setw(4) << atom1 << '_' << xyz[axis] << '_' << (sign > 0 ? 'p' : 'm');
    return o.str();
}

} // namespace fdvib
