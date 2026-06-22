#include "fdvib.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

namespace fdvib {

namespace {

bool has_pw_output(const fs::path &workdir) {
    if (!fs::exists(workdir)) return false;
    for (const auto &e : fs::recursive_directory_iterator(workdir))
        if (e.is_regular_file() && e.path().extension() == ".out" &&
            e.path().filename().string().rfind("disp_", 0) == 0) return true;
    return false;
}

void snapshot_or_check(const fs::path &src, const fs::path &dst,
                       const fs::path &workdir, bool force) {
    const auto wanted = read_text(src);
    if (fs::exists(dst)) {
        if (read_text(dst) != wanted) {
            if (!force || has_pw_output(workdir))
                throw std::runtime_error("Existing snapshot differs; --force is allowed only before any pw.x output: " + dst.string());
            write_text(dst, wanted);
        }
    } else write_text(dst, wanted);
}

} // namespace

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
    std::string line; std::smatch m; bool started = false; int count = 0;
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
            write_forces(forces, parse_forces(output, q.nat), output);
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

} // namespace fdvib
