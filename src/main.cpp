#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#ifndef FDVIB_VERSION
#define FDVIB_VERSION "unknown"
#endif

namespace fs = std::filesystem;

constexpr double BOHR_TO_ANG = 0.529177210903;
constexpr double RY_TO_EV = 13.605693122994;
// QE Rydberg atomic units use m_e = 1/2, hence AMU_RY = AMU_AU / 2.
constexpr double AMU_RY = 911.4442431045;
constexpr double CM_TO_EV = 1.239841984e-4;
constexpr double KB_EV = 8.617333262145e-5;
constexpr double KB_SI = 1.380649e-23;
constexpr double H_SI = 6.62607015e-34;
constexpr double AMU_KG = 1.66053906660e-27;
constexpr double BAR_PA = 1.0e5;
constexpr double PI = 3.141592653589793238462643383279502884;

using Vec3 = std::array<double, 3>;

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

struct Config {
    std::map<std::string, std::string> v;

    static Config load(const fs::path &p) {
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

    bool has(const std::string &k) const { return v.count(lower(k)) != 0; }
    std::string get(const std::string &k, const std::string &d = {}) const {
        const auto it = v.find(lower(k));
        return it == v.end() ? d : unquote(it->second);
    }
    double real(const std::string &k, double d) const {
        return has(k) ? number(v.at(lower(k))) : d;
    }
    int integer(const std::string &k, int d) const {
        if (!has(k)) return d;
        const double x = number(v.at(lower(k)));
        const double rounded = std::round(x);
        if (std::abs(x - rounded) > 1.0e-10)
            throw std::runtime_error("Expected integer for " + k + ": " + v.at(lower(k)));
        return static_cast<int>(rounded);
    }

    void require_only(const std::set<std::string> &allowed, const std::string &context) const {
        for (const auto &[key, value] : v) {
            (void)value;
            if (!allowed.count(key))
                throw std::runtime_error("Unknown parameter in " + context + ": " + key);
        }
    }
};

std::vector<int> integer_list(std::string s) {
    for (char &c : s) if (c == ',' || c == ';') c = ' ';
    std::istringstream in(s);
    std::vector<int> out;
    int x;
    while (in >> x) out.push_back(x);
    return out;
}

struct Species { std::string symbol, pseudo; double mass{}; };
struct Atom { std::string symbol; Vec3 r{}; std::vector<std::string> extra; int type{}; };

struct QEInput {
    std::string text, clean_text;
    std::vector<std::string> lines;
    int nat{}, ntyp{}, pos_header{}, pos_start{}, cell_header{};
    std::vector<Species> species;
    std::vector<Atom> atoms;
    std::array<Vec3, 3> cell{};
};

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
    if (std::stoi(mib[1]) != 0) throw std::runtime_error("First version requires ibrav=0");
    if (!std::regex_search(q.clean_text, std::regex(R"(\bcalculation\s*=\s*['\"]scf['\"])", std::regex::icase)))
        throw std::runtime_error("scf.in must contain calculation='scf'");
    if (!std::regex_search(q.clean_text, std::regex(R"(\btprnfor\s*=\s*\.true\.)", std::regex::icase)))
        throw std::runtime_error("scf.in must contain tprnfor=.true.");

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

std::string displaced_input(const QEInput &q, int atom, int axis, double shift,
                            const std::string &outdir) {
    auto lines = q.lines;
    auto a = q.atoms.at(atom);
    a.r[axis] += shift;
    lines[q.pos_start + atom] = format_position(a);
    bool found = false;
    std::regex ore(R"(^\s*outdir\s*=)", std::regex::icase);
    for (auto &l : lines) {
        if (std::regex_search(l, ore)) {
            const auto indent = l.substr(0, l.find_first_not_of(" \t"));
            l = indent + "outdir = '" + outdir + "'\n";
            found = true;
        }
    }
    if (!found) throw std::runtime_error("Cannot find outdir in scf.in");
    return std::accumulate(lines.begin(), lines.end(), std::string{});
}

struct Settings {
    fs::path config_path, root, scf_input, workdir;
    std::string system_type, pw_command, mpi_command, output_prefix;
    std::vector<int> selected;
    double displacement{};
    int multiplicity{};
    bool multiplicity_explicit{}, selected_all{};
};

Settings settings(const fs::path &config_path, const fs::path &root_override = {}) {
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

bool has_pw_output(const fs::path &workdir) {
    if (!fs::exists(workdir)) return false;
    for (const auto &e : fs::recursive_directory_iterator(workdir))
        if (e.is_regular_file() && e.path().extension() == ".out" &&
            e.path().filename().string().rfind("disp_", 0) == 0) return true;
    return false;
}

void snapshot_or_check(const fs::path &src, const fs::path &dst, const fs::path &workdir, bool force) {
    const auto wanted = read_text(src);
    if (fs::exists(dst)) {
        if (read_text(dst) != wanted) {
            if (!force || has_pw_output(workdir))
                throw std::runtime_error("Existing snapshot differs; --force is allowed only before any pw.x output: " + dst.string());
            write_text(dst, wanted);
        }
    } else write_text(dst, wanted);
}

void prepare(const Settings &s, bool force) {
    const auto q = parse_qe_input(s.scf_input);
    if (s.system_type == "gas") {
        if (!s.multiplicity_explicit || s.multiplicity < 1)
            throw std::runtime_error("gas requires an explicit positive multiplicity in fdvib.in");
        if (!s.selected_all)
            throw std::runtime_error("gas requires selected_atoms='all'");
        std::smatch m;
        int nspin = 1;
        if (std::regex_search(q.clean_text, m, std::regex(R"(\bnspin\s*=\s*(\d+))", std::regex::icase)))
            nspin = std::stoi(m[1]);
        double magnetization = 0.0;
        const bool has_magnetization = std::regex_search(
            q.clean_text, m, std::regex(R"(\btot_magnetization\s*=\s*([-+0-9.EeDd]+))", std::regex::icase));
        if (has_magnetization) magnetization = number(m[1]);
        if (s.multiplicity == 1) {
            if (nspin != 1 && !(nspin == 2 && has_magnetization && std::abs(magnetization) < 1e-8))
                throw std::runtime_error("Gas singlet requires nspin=1, or nspin=2 with tot_magnetization=0");
        } else if (nspin != 2 || !has_magnetization ||
                   std::abs(magnetization - (s.multiplicity - 1)) > 1e-8) {
            throw std::runtime_error("Gas multiplicity requires nspin=2 and tot_magnetization=multiplicity-1");
        }
    }
    auto selected = s.selected;
    if (s.system_type == "gas") {
        if (!selected.empty() && static_cast<int>(selected.size()) != q.nat)
            throw std::runtime_error("gas requires selected_atoms='all'");
        selected.resize(q.nat); std::iota(selected.begin(), selected.end(), 1);
    }
    if (selected.empty()) throw std::runtime_error("selected_atoms is empty");
    std::set<int> unique;
    for (int a : selected) {
        if (a < 1 || a > q.nat || !unique.insert(a).second) throw std::runtime_error("Invalid/duplicate selected atom");
    }
    fs::create_directories(s.workdir);
    snapshot_or_check(s.scf_input, s.workdir / "scf.in.reference", s.workdir, force);
    snapshot_or_check(s.config_path, s.workdir / "fdvib.in.reference", s.workdir, force);

    std::ostringstream jobs;
    jobs << "# job atom axis sign displacement_A directory\n";
    for (int atom1 : selected) for (int axis = 0; axis < 3; ++axis) for (int sign : {1, -1}) {
        const auto name = job_name(atom1, axis, sign);
        const auto dir = s.workdir / name;
        const auto input = dir / (name + ".in");
        const auto output = dir / (name + ".out");
        const auto relout = fs::relative(dir / "out", s.root).generic_string();
        const auto wanted = displaced_input(q, atom1 - 1, axis, sign * s.displacement, "./" + relout);
        fs::create_directories(dir);
        if (fs::exists(input)) {
            if (read_text(input) != wanted) {
                if (fs::exists(output)) throw std::runtime_error("Refuse to modify a completed/started task: " + name);
                if (!force) throw std::runtime_error("Existing input differs (use --force only before run): " + name);
                write_text(input, wanted);
            }
        } else write_text(input, wanted);
        static const char *xyz[] = {"x", "y", "z"};
        jobs << name << ' ' << atom1 << ' ' << xyz[axis] << ' ' << (sign > 0 ? '+' : '-') << ' '
             << std::fixed << std::setprecision(8) << sign * s.displacement << ' ' << name << '\n';
    }
    const auto jp = s.workdir / "jobs.list";
    if (fs::exists(jp) && read_text(jp) != jobs.str() && !force)
        throw std::runtime_error("Existing jobs.list differs");
    if (!fs::exists(jp) || read_text(jp) != jobs.str()) write_text(jp, jobs.str());
    std::cout << "Prepared " << selected.size() * 6 << " displacement jobs in " << s.workdir << '\n';
}

std::string shell_quote(const std::string &s) {
    std::string q = "'";
    for (char c : s) q += (c == '\'') ? "'\\''" : std::string(1, c);
    return q + "'";
}

int shell_run(const std::string &cmd, const fs::path &cwd, const fs::path &stdout_path) {
    const pid_t pid = fork();
    if (pid < 0) throw std::runtime_error("fork failed");
    if (pid == 0) {
        if (chdir(cwd.c_str()) != 0) _exit(126);
        FILE *f = std::fopen(stdout_path.c_str(), "w");
        if (!f) _exit(126);
        dup2(fileno(f), STDOUT_FILENO); dup2(fileno(f), STDERR_FILENO); std::fclose(f);
        execl("/bin/sh", "sh", "-c", cmd.c_str(), static_cast<char *>(nullptr));
        _exit(127);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return 128 + (WIFSIGNALED(status) ? WTERMSIG(status) : 0);
}

std::vector<Vec3> parse_forces(const fs::path &output, int nat) {
    const auto text = read_text(output);
    if (text.find("JOB DONE") == std::string::npos) throw std::runtime_error("JOB DONE not found: " + output.string());
    if (std::regex_search(text, std::regex(R"(convergence\s+NOT\s+achieved)", std::regex::icase)))
        throw std::runtime_error("SCF convergence was not achieved: " + output.string());
    if (std::regex_search(text, std::regex(R"(Error\s+in\s+routine)", std::regex::icase)))
        throw std::runtime_error("QE reported Error in routine: " + output.string());
    const auto pos = text.rfind("Forces acting on atoms");
    if (pos == std::string::npos) throw std::runtime_error("Force block not found: " + output.string());
    std::istringstream in(text.substr(pos));
    std::regex re(R"(^\s*atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+))", std::regex::icase);
    std::vector<Vec3> f(nat);
    std::vector<bool> got(nat, false);
    std::string line; std::smatch m;
    while (std::getline(in, line)) if (std::regex_search(line, m, re)) {
        const int i = std::stoi(m[1]) - 1;
        if (i >= 0 && i < nat) {
            f[i] = {number(m[2]), number(m[3]), number(m[4])}; got[i] = true;
        }
    }
    if (std::count(got.begin(), got.end(), true) != nat) throw std::runtime_error("Incomplete force block: " + output.string());
    return f;
}

void write_forces(const fs::path &p, const std::vector<Vec3> &f, const fs::path &source) {
    std::ostringstream o;
    o << "# source: " << source.filename().string() << "\n# nat: " << f.size() << "\n# units: Ry/Bohr\n# atom fx fy fz\n";
    o << std::scientific << std::setprecision(15);
    for (std::size_t i = 0; i < f.size(); ++i) o << i + 1 << ' ' << f[i][0] << ' ' << f[i][1] << ' ' << f[i][2] << '\n';
    write_text(p, o.str());
}

std::vector<Vec3> read_forces(const fs::path &p, int nat) {
    std::istringstream in(read_text(p));
    std::vector<Vec3> f(nat); std::vector<bool> got(nat, false);
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line); if (line.empty() || line.front() == '#') continue;
        std::istringstream ls(line); int i; Vec3 v;
        if (!(ls >> i >> v[0] >> v[1] >> v[2]) || i < 1 || i > nat) throw std::runtime_error("Bad forces.dat: " + p.string());
        f[i - 1] = v; got[i - 1] = true;
    }
    if (std::count(got.begin(), got.end(), true) != nat) throw std::runtime_error("Incomplete forces.dat: " + p.string());
    return f;
}

void run_jobs(const Settings &s) {
    const auto q = parse_qe_input(s.workdir / "scf.in.reference");
    std::istringstream jobs(read_text(s.workdir / "jobs.list"));
    std::string line;
    int done = 0, skipped = 0;
    while (std::getline(jobs, line)) {
        line = trim(line); if (line.empty() || line.front() == '#') continue;
        std::istringstream ls(line); std::string id, axis, sign, rel; int atom; double disp;
        ls >> id >> atom >> axis >> sign >> disp >> rel;
        const auto dir = s.workdir / rel, input = dir / (id + ".in"), output = dir / (id + ".out"), forces = dir / "forces.dat";
        if (fs::exists(output)) {
            if (!fs::exists(forces)) write_forces(forces, parse_forces(output, q.nat), output);
            ++skipped; continue;
        }
        const std::string cmd = (s.mpi_command.empty() ? "" : s.mpi_command + " ") + s.pw_command + " -in " + shell_quote(fs::relative(input, s.root).string());
        std::cout << "Running " << id << std::endl;
        const int rc = shell_run(cmd, s.root, fs::relative(output, s.root));
        std::ostringstream status; status << "exit_code " << rc << '\n';
        write_text(dir / "status", status.str());
        if (rc != 0) throw std::runtime_error(id + " failed with exit code " + std::to_string(rc));
        write_forces(forces, parse_forces(output, q.nat), output);
        ++done;
    }
    std::cout << "Completed " << done << ", preserved " << skipped << " existing jobs\n";
}

double norm(const Vec3 &v) { return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }

void write_dynG(const fs::path &p, const QEInput &q, const std::vector<double> &h) {
    const double alat_ang = norm(q.cell[0]);
    const double alat_bohr = alat_ang / BOHR_TO_ANG;
    std::ostringstream o;
    o << "Dynamical matrix file\nfdvib finite-difference local/Gamma matrix\n";
    o << std::setw(3) << q.ntyp << std::setw(5) << q.nat << std::setw(4) << 0;
    o << std::fixed << std::setprecision(7) << std::setw(12) << alat_bohr;
    for (int i = 0; i < 5; ++i) o << std::setw(12) << 0.0;
    o << "\nBasis vectors\n" << std::setprecision(9);
    for (const auto &v : q.cell) o << "  " << std::setw(15) << v[0]/alat_ang << std::setw(15) << v[1]/alat_ang << std::setw(15) << v[2]/alat_ang << '\n';
    o << std::setprecision(12);
    for (int i = 0; i < q.ntyp; ++i)
        o << std::setw(12) << i+1 << "  '" << std::left << std::setw(7) << q.species[i].symbol << std::right << "' " << std::setw(22) << AMU_RY*q.species[i].mass << '\n';
    o << std::setprecision(10);
    for (int i = 0; i < q.nat; ++i)
        o << std::setw(5) << i+1 << std::setw(5) << q.atoms[i].type << std::setw(18) << q.atoms[i].r[0]/alat_ang
          << std::setw(18) << q.atoms[i].r[1]/alat_ang << std::setw(18) << q.atoms[i].r[2]/alat_ang << '\n';
    o << "\n     Dynamical  Matrix in cartesian axes\n\n     q = (    0.000000000   0.000000000   0.000000000 ) \n\n";
    o << std::fixed << std::setprecision(10);
    const int n3 = 3*q.nat;
    for (int a = 0; a < q.nat; ++a) for (int b = 0; b < q.nat; ++b) {
        o << std::setw(5) << a+1 << std::setw(5) << b+1 << '\n';
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) o << std::setw(14) << h[(3*a+i)*n3 + 3*b+j] << "   " << std::setw(12) << 0.0 << "  ";
            o << '\n';
        }
    }
    write_text(p, o.str());
}

struct Mode { double freq{}; std::vector<Vec3> displacement; };

std::vector<Mode> parse_modes(const fs::path &p, int nat) {
    std::ifstream in(p);
    if (!in) throw std::runtime_error("Cannot read " + p.string());
    std::regex fr(R"(freq\s*\(\s*\d+\s*\)\s*=\s*[-+0-9.EeDd]+\s*\[THz\]\s*=\s*([-+0-9.EeDd]+)\s*\[cm-1\])", std::regex::icase);
    std::regex pair(R"(\(\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*\))");
    std::vector<Mode> out; std::string line; std::smatch m;
    while (std::getline(in, line)) if (std::regex_search(line, m, fr)) {
        Mode mode; mode.freq = number(m[1]); mode.displacement.reserve(nat);
        for (int i = 0; i < nat; ++i) {
            if (!std::getline(in, line) || !std::regex_search(line, m, pair)) throw std::runtime_error("Bad eigenvector block in " + p.string());
            mode.displacement.push_back({number(m[1]), number(m[3]), number(m[5])});
        }
        out.push_back(std::move(mode));
    }
    if (static_cast<int>(out.size()) != 3 * nat)
        throw std::runtime_error("Expected " + std::to_string(3 * nat) + " modes, found " +
                                 std::to_string(out.size()) + " in " + p.string());
    return out;
}

struct DynGeometry;

void write_compact_molden(const fs::path &p, const DynGeometry &g,
                          const std::vector<Mode> &modes, double zero_tol);

void analyze(const Settings &s) {
    const auto q = parse_qe_input(s.workdir / "scf.in.reference");
    auto selected = s.selected;
    if (s.system_type == "gas") { selected.resize(q.nat); std::iota(selected.begin(), selected.end(), 1); }
    const int n3 = 3*q.nat;
    std::vector<double> h(n3*n3, 0.0);
    const double delta_bohr = s.displacement / BOHR_TO_ANG;
    for (int atom1 : selected) for (int axis = 0; axis < 3; ++axis) {
        const auto p = read_forces(s.workdir/job_name(atom1,axis,1)/"forces.dat", q.nat);
        const auto m = read_forces(s.workdir/job_name(atom1,axis,-1)/"forces.dat", q.nat);
        const int col = 3*(atom1-1)+axis;
        for (int atomj1 : selected) for (int beta=0; beta<3; ++beta) {
            const int row=3*(atomj1-1)+beta;
            h[row*n3+col]=-(p[atomj1-1][beta]-m[atomj1-1][beta])/(2*delta_bohr);
        }
    }
    double asym=0;
    for (int i=0;i<n3;++i) for(int j=0;j<i;++j) {
        asym=std::max(asym,std::abs(h[i*n3+j]-h[j*n3+i]));
        const double x=0.5*(h[i*n3+j]+h[j*n3+i]); h[i*n3+j]=h[j*n3+i]=x;
    }
    const auto results=s.workdir/"results"; fs::create_directories(results);
    const auto dyn=results/(s.output_prefix+".dynG");
    for (const auto &p : {dyn, results/"dynmat.in"})
        if (fs::exists(p)) throw std::runtime_error("Refuse to overwrite " + p.string());
    const auto result_config = results / "fdvib.in.reference";
    if (fs::exists(result_config)) {
        if (read_text(result_config) != read_text(s.config_path))
            throw std::runtime_error("Result configuration snapshot differs: " + result_config.string());
    } else {
        fs::copy_file(s.config_path, result_config);
    }
    write_dynG(dyn,q,h);
    std::ostringstream di;
    di << "&INPUT\n  fildyn='" << dyn.filename().string() << "',\n  filout='" << s.output_prefix << ".freq.out',\n"
       << "  filmol=' ',\n  asr='" << (s.system_type=="gas"?"zero-dim":"no") << "',\n"
       << "  remove_interaction_blocks=" << (s.system_type=="gas"?".false.":".true.") << ",\n/\n";
    write_text(results/"dynmat.in",di.str());
    const auto root_thermo=s.root/"thermo.in";
    if(fs::exists(root_thermo)&&!fs::exists(results/"thermo.in")) fs::copy_file(root_thermo,results/"thermo.in");
    std::cout << "Hessian max asymmetry: " << std::scientific << asym << " Ry/Bohr^2\n"
              << "Wrote " << dyn << " and " << (results/"dynmat.in") << "\n"
              << "Next, run in the result directory:\n  dynmat.x -in dynmat.in > dynmat.out\n";
}

struct DynGeometry { std::vector<double> masses; std::vector<Vec3> r_bohr; std::vector<std::string> symbols; };

DynGeometry read_dyn_geometry(const fs::path &p) {
    std::ifstream in(p); if(!in) throw std::runtime_error("Cannot read "+p.string());
    std::string line; std::getline(in,line); std::getline(in,line); std::getline(in,line);
    std::istringstream hs(line); int ntyp,nat,ibrav; double alat; hs>>ntyp>>nat>>ibrav>>alat;
    if(!hs||ibrav!=0) throw std::runtime_error("Unsupported dynG header");
    std::getline(in,line); for(int i=0;i<3;++i) std::getline(in,line);
    std::vector<double> tmass(ntyp); std::vector<std::string> tsymbol(ntyp);
    const std::regex quoted(R"('([^']*)')");
    for(int i=0;i<ntyp;++i){
        std::getline(in,line);
        std::smatch match;
        if (!std::regex_search(line, match, quoted)) throw std::runtime_error("Bad species symbol in dynG");
        tsymbol[i]=trim(match[1]);
        std::istringstream ls(line); std::string token, last;
        while(ls>>token) last=token;
        if(last.empty()) throw std::runtime_error("Bad species mass in dynG");
        tmass[i]=number(last)/AMU_RY;
    }
    DynGeometry g; g.masses.resize(nat); g.r_bohr.resize(nat); g.symbols.resize(nat);
    for(int i=0;i<nat;++i){ std::getline(in,line); std::istringstream ls(line); int n,t; Vec3 x; ls>>n>>t>>x[0]>>x[1]>>x[2]; g.masses[i]=tmass.at(t-1);g.symbols[i]=tsymbol.at(t-1);for(int k=0;k<3;++k)g.r_bohr[i][k]=x[k]*alat; }
    return g;
}

void write_compact_molden(const fs::path &p, const DynGeometry &g,
                          const std::vector<Mode> &modes, double zero_tol) {
    std::vector<const Mode*> keep;
    for (const auto &m : modes) {
        double n = 0;
        for (const auto &d : m.displacement) for (double x : d) n += x*x;
        if (std::abs(m.freq) >= zero_tol && n > 1e-16) keep.push_back(&m);
    }
    std::ostringstream o;
    o << "[Molden Format]\n[FREQ]\n" << std::fixed << std::setprecision(8);
    for (auto m : keep) o << m->freq << '\n';
    o << "[FR-COORD]\n";
    for (std::size_t i=0;i<g.symbols.size();++i) o << std::setw(6) << g.symbols[i] << ' ' << std::setw(15) << g.r_bohr[i][0] << ' '
                                    << std::setw(15) << g.r_bohr[i][1] << ' ' << std::setw(15) << g.r_bohr[i][2] << '\n';
    o << "[FR-NORM-COORD]\n" << std::setprecision(10);
    for (std::size_t i = 0; i < keep.size(); ++i) {
        o << " vibration " << i+1 << '\n';
        for (const auto &d : keep[i]->displacement) o << std::setw(15) << d[0] << ' ' << std::setw(15) << d[1] << ' ' << std::setw(15) << d[2] << '\n';
    }
    write_text(p, o.str());
    std::cout << "Wrote " << keep.size() << " compact Molden modes\n";
}

void modes(const fs::path &results) {
    if (!fs::is_directory(results)) throw std::runtime_error("Not a result directory: " + results.string());
    std::vector<fs::path> dyns, freqs;
    for(const auto&e:fs::directory_iterator(results)){const auto n=e.path().filename().string();if(n.size()>=5&&n.substr(n.size()-5)==".dynG")dyns.push_back(e.path());if(n.size()>=9&&n.substr(n.size()-9)==".freq.out")freqs.push_back(e.path());}
    if(dyns.size()!=1||freqs.size()!=1) throw std::runtime_error("Result directory must contain exactly one .dynG and one .freq.out");
    const auto &dyn=dyns.front(), &freq=freqs.front();
    const auto g=read_dyn_geometry(dyn);
    const auto parsed=parse_modes(freq,static_cast<int>(g.masses.size()));
    const auto mold=results/(dyn.stem().string()+".mold");
    if(fs::exists(mold)) throw std::runtime_error("Refuse to overwrite "+mold.string());
    write_compact_molden(mold,g,parsed,1.0e-6);
}

std::array<double,3> symmetric_eigenvalues(Vec3 diag, Vec3 off) {
    // Stable analytic eigenvalues of a real symmetric 3x3 matrix.
    const double p1=off[0]*off[0]+off[1]*off[1]+off[2]*off[2];
    if(p1==0){ std::array<double,3>a{diag[0],diag[1],diag[2]}; std::sort(a.begin(),a.end()); return a; }
    const double q=(diag[0]+diag[1]+diag[2])/3;
    const double p2=(diag[0]-q)*(diag[0]-q)+(diag[1]-q)*(diag[1]-q)+(diag[2]-q)*(diag[2]-q)+2*p1;
    const double p=std::sqrt(p2/6);
    const double a00=(diag[0]-q)/p,a11=(diag[1]-q)/p,a22=(diag[2]-q)/p;
    const double a01=off[0]/p,a02=off[1]/p,a12=off[2]/p;
    const double det=a00*a11*a22+2*a01*a02*a12-a00*a12*a12-a11*a02*a02-a22*a01*a01;
    const double r=std::clamp(det/2.0,-1.0,1.0), phi=std::acos(r)/3;
    const double e3=q+2*p*std::cos(phi), e1=q+2*p*std::cos(phi+2*PI/3), e2=3*q-e1-e3;
    std::array<double,3>a{e1,e2,e3}; std::sort(a.begin(),a.end()); return a;
}

struct VibThermo { double zpe{},u{},s{},f{}; int imag{},zero{},floored{},used{}; };

VibThermo vibration_thermo(const std::vector<Mode>& modes,double T,const std::string& model,double floor,double zero_tol) {
    VibThermo v;
    for(const auto&m:modes){ double nu=m.freq;
        if(std::abs(nu)<zero_tol){++v.zero;continue;} if(nu<0){++v.imag;continue;}
        if(model=="frequency_floor"&&nu<floor){nu=floor;++v.floored;} ++v.used;
        const double e=nu*CM_TO_EV,x=e/(KB_EV*T);
        v.zpe+=0.5*e; v.u+=0.5*e+e/std::expm1(x);
        v.s+=KB_EV*(x/std::expm1(x)-std::log1p(-std::exp(-x)));
    }
    v.f=v.u-T*v.s; return v;
}

void thermo(const fs::path &results) {
    if (!fs::is_directory(results)) throw std::runtime_error("Not a result directory: " + results.string());
    const auto c=Config::load(results/"thermo.in");
    c.require_only({"model", "temperature_k", "pressure_bar", "symmetry_number",
                    "electronic_degeneracy", "rotor_type", "low_frequency_model",
                    "frequency_floor_cm1", "zero_tolerance_cm1"}, "thermo.in");
    const std::string model=lower(c.get("model"));
    if(model!="gas_rrho"&&model!="local_harmonic") throw std::runtime_error("THERMO model must be gas_rrho or local_harmonic");
    if(model=="local_harmonic" &&
       (c.has("pressure_bar") || c.has("symmetry_number") ||
        c.has("electronic_degeneracy") || c.has("rotor_type")))
        throw std::runtime_error("local_harmonic must not contain gas-only thermochemistry parameters");
    if (fs::exists(results/"dynmat.in")) {
        const auto d=Config::load(results/"dynmat.in");
        const auto asr=lower(d.get("asr","no"));
        const auto remove=lower(d.get("remove_interaction_blocks",".false."));
        const bool removes_blocks=(remove==".true."||remove=="true"||remove=="t");
        if(model=="local_harmonic"&&!removes_blocks)
            throw std::runtime_error("local_harmonic requires remove_interaction_blocks=.true. in dynmat.in");
        if(model=="gas_rrho"&&(removes_blocks||asr!="zero-dim"))
            throw std::runtime_error("gas_rrho requires asr='zero-dim' and remove_interaction_blocks=.false. in dynmat.in");
    }
    const double T=c.real("temperature_k",-1); if(T<=0) throw std::runtime_error("temperature_k must be > 0");
    const std::string low=lower(c.get("low_frequency_model","harmonic"));
    if(low!="harmonic"&&low!="frequency_floor") throw std::runtime_error("Bad low_frequency_model");
    const double floor=c.real("frequency_floor_cm1",50),zt=c.real("zero_tolerance_cm1",1);
    if (zt < 0.0) throw std::runtime_error("zero_tolerance_cm1 must be non-negative");
    if (low=="frequency_floor" && floor <= 0.0) throw std::runtime_error("frequency_floor_cm1 must be positive");
    std::vector<fs::path> dyns, freqs;
    for(const auto&e:fs::directory_iterator(results)){const auto n=e.path().filename().string();if(n.size()>=5&&n.substr(n.size()-5)==".dynG")dyns.push_back(e.path());if(n.size()>=9&&n.substr(n.size()-9)==".freq.out")freqs.push_back(e.path());}
    if(dyns.size()!=1||freqs.size()!=1) throw std::runtime_error("Result directory must contain exactly one .dynG and one .freq.out");
    const auto &dyn=dyns.front(), &freq=freqs.front();
    const auto modes=parse_modes(freq,read_dyn_geometry(dyn).masses.size());
    const auto vib=vibration_thermo(modes,T,low,floor,zt);
    double hcorr=vib.u,stotal=vib.s,gcorr=vib.f,strans=0,srot=0,selec=0,htrans=0,urot=0;
    if(model=="gas_rrho"){
        if(!c.has("pressure_bar")||!c.has("symmetry_number")) throw std::runtime_error("gas_rrho requires explicit pressure_bar and symmetry_number");
        const auto dataset=Config::load(results/"fdvib.in.reference");
        if(lower(dataset.get("system_type"))!="gas") throw std::runtime_error("gas_rrho requires a gas fdvib.in.reference");
        const double pbar=c.real("pressure_bar",-1); const int sigma=c.integer("symmetry_number",0),mult=dataset.integer("multiplicity",0);
        if(pbar<=0||sigma<=0||mult<=0) throw std::runtime_error("Gas pressure, symmetry, and multiplicity must be positive");
        const auto g=read_dyn_geometry(dyn); const double mtot=std::accumulate(g.masses.begin(),g.masses.end(),0.0);
        Vec3 com{0,0,0}; for(size_t i=0;i<g.masses.size();++i)for(int k=0;k<3;++k)com[k]+=g.masses[i]*g.r_bohr[i][k]/mtot;
        Vec3 d{0,0,0},off{0,0,0};
        for(size_t i=0;i<g.masses.size();++i){const double m=g.masses[i]*AMU_KG;double x[3];for(int k=0;k<3;++k)x[k]=(g.r_bohr[i][k]-com[k])*BOHR_TO_ANG*1e-10;
            d[0]+=m*(x[1]*x[1]+x[2]*x[2]);d[1]+=m*(x[0]*x[0]+x[2]*x[2]);d[2]+=m*(x[0]*x[0]+x[1]*x[1]);off[0]-=m*x[0]*x[1];off[1]-=m*x[0]*x[2];off[2]-=m*x[1]*x[2];}
        const auto I=symmetric_eigenvalues(d,off); std::string type=lower(c.get("rotor_type","auto"));
        if(type=="auto") type=(g.masses.size()==1?"atom":(I[0]/std::max(I[2],1e-300)<1e-6?"linear":"nonlinear"));
        const double qtrans=std::pow(2*PI*(mtot*AMU_KG)*KB_SI*T/(H_SI*H_SI),1.5)*(KB_SI*T/(pbar*BAR_PA));
        strans=KB_EV*(std::log(qtrans)+2.5);htrans=2.5*KB_EV*T;
        if(type=="linear") {
            const double moment=std::max(I[1],I[2]); if(!(moment>0.0)) throw std::runtime_error("Invalid linear-molecule moment of inertia");
            const double qr=8*PI*PI*moment*KB_SI*T/(sigma*H_SI*H_SI);srot=KB_EV*(std::log(qr)+1);urot=KB_EV*T;
        }
        else if(type=="nonlinear") {
            if(!(I[0]>0.0&&I[1]>0.0&&I[2]>0.0)) throw std::runtime_error("Invalid nonlinear-molecule moments of inertia");
            const double qr=std::sqrt(PI)/sigma*std::pow(8*PI*PI*KB_SI*T/(H_SI*H_SI),1.5)*std::sqrt(I[0]*I[1]*I[2]);srot=KB_EV*(std::log(qr)+1.5);urot=1.5*KB_EV*T;
        }
        else if(type!="atom") throw std::runtime_error("rotor_type must be auto, atom, linear, nonlinear");
        int degeneracy=mult; const auto ed=lower(c.get("electronic_degeneracy","auto"));
        if(ed!="auto") {
            const double value=number(ed), rounded=std::round(value);
            if(std::abs(value-rounded)>1.0e-10) throw std::runtime_error("electronic_degeneracy must be an integer or auto");
            degeneracy=static_cast<int>(rounded);
        }
        if(degeneracy<=0)throw std::runtime_error("electronic_degeneracy must be positive or auto");
        selec=KB_EV*std::log(static_cast<double>(degeneracy)); hcorr=vib.u+htrans+urot; stotal=vib.s+strans+srot+selec; gcorr=hcorr-T*stotal;
    }
    std::ostringstream o;o<<"# FDVIB thermochemistry\n# model: "<<model<<"\n# low_frequency_model: "<<low<<"\n";
    if(low=="frequency_floor")o<<"# frequency_floor_cm1: "<<floor<<"\n";
    o<<"# imaginary_modes_excluded: "<<vib.imag<<"\n# zero_modes_excluded: "<<vib.zero<<"\n# positive_modes_used: "<<vib.used<<"\n# modes_floored: "<<vib.floored<<"\n";
    o<<"# units: T=K energies=eV entropy=eV/K\n"<<std::fixed;
    if(model=="local_harmonic") {
        o<<"# "<<std::setw(11)<<"T/K"<<std::setw(16)<<"ZPE/eV"<<std::setw(16)<<"U_vib/eV"
         <<std::setw(19)<<"S_vib/eV_K"<<std::setw(16)<<"TS_vib/eV"<<std::setw(16)<<"F_vib/eV"<<'\n';
        o<<"  "<<std::setw(11)<<std::setprecision(3)<<T
         <<std::setw(16)<<std::setprecision(10)<<vib.zpe<<std::setw(16)<<vib.u
         <<std::setw(19)<<std::setprecision(12)<<vib.s
         <<std::setw(16)<<std::setprecision(10)<<T*vib.s<<std::setw(16)<<vib.f<<'\n';
    } else {
        o<<"# "<<std::setw(11)<<"T/K"<<std::setw(15)<<"ZPE/eV"<<std::setw(15)<<"U_vib/eV"
         <<std::setw(15)<<"H_trans/eV"<<std::setw(15)<<"U_rot/eV"
         <<std::setw(18)<<"S_trans/eV_K"<<std::setw(18)<<"S_rot/eV_K"
         <<std::setw(18)<<"S_vib/eV_K"<<std::setw(18)<<"S_elec/eV_K"
         <<std::setw(15)<<"H_corr/eV"<<std::setw(15)<<"G_corr/eV"<<'\n';
        o<<"  "<<std::setw(11)<<std::setprecision(3)<<T
         <<std::setw(15)<<std::setprecision(10)<<vib.zpe<<std::setw(15)<<vib.u
         <<std::setw(15)<<htrans<<std::setw(15)<<urot
         <<std::setw(18)<<std::setprecision(12)<<strans<<std::setw(18)<<srot
         <<std::setw(18)<<vib.s<<std::setw(18)<<selec
         <<std::setw(15)<<std::setprecision(10)<<hcorr<<std::setw(15)<<gcorr<<'\n';
    }
    write_text(results/"thermo.dat",o.str());
    std::cout<<o.str()<<"Written "<<(results/"thermo.dat")<<'\n';
}

void usage() {
    std::cerr << "Usage:\n  fdvib prepare fdvib.in [--force]\n  fdvib run fdvib.in\n  fdvib analyze fdvib.in\n  fdvib modes RESULTS_DIR\n  fdvib thermo RESULTS_DIR\n";
}

int main(int argc,char**argv){
    try{
        if(argc==2&&std::string(argv[1])=="--help"){usage();return 0;}
        if(argc==2&&std::string(argv[1])=="--version"){std::cout<<"fdvib "<<FDVIB_VERSION<<'\n';return 0;}
        if(argc<3){usage();return 2;} const std::string cmd=argv[1];
        if(cmd=="prepare") prepare(settings(argv[2]),argc==4&&std::string(argv[3])=="--force");
        else if(cmd=="run") run_jobs(settings(argv[2]));
        else if(cmd=="analyze") {
            const auto config_path = fs::absolute(argv[2]);
            const auto root = config_path.parent_path();
            const auto current = Config::load(config_path);
            const auto workdir = root / current.get("workdir", "fdvib");
            analyze(settings(workdir / "fdvib.in.reference", root));
        }
        else if(cmd=="modes") modes(fs::absolute(argv[2]));
        else if(cmd=="thermo") thermo(fs::absolute(argv[2]));
        else {usage();return 2;}
        return 0;
    }catch(const std::exception&e){std::cerr<<"fdvib: error: "<<e.what()<<'\n';return 1;}
}
