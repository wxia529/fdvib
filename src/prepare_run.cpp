#include "fdvib.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <cstdint>
#include <stdexcept>
#include <sys/wait.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>

namespace fdvib {

namespace {

class CalculationLock {
  public:
    explicit CalculationLock(const fs::path &workdir) {
        path_ = workdir.string() + ".lock";
        fd_ = open(path_.c_str(), O_CREAT | O_RDWR, 0644);
        if (fd_ < 0) throw std::runtime_error("Cannot open calculation lock: " + path_.string());
        if (flock(fd_, LOCK_EX | LOCK_NB) != 0) {
            close(fd_);
            fd_ = -1;
            throw std::runtime_error("Calculation is already running: " + workdir.string());
        }
        const auto metadata = "pid=" + std::to_string(static_cast<long long>(getpid())) + "\n";
        if (ftruncate(fd_, 0) == 0) (void)write(fd_, metadata.data(), metadata.size());
    }
    ~CalculationLock() {
        if (fd_ >= 0) {
            (void)flock(fd_, LOCK_UN);
            close(fd_);
        }
    }
    CalculationLock(const CalculationLock &) = delete;
    CalculationLock &operator=(const CalculationLock &) = delete;
  private:
    int fd_{-1};
    fs::path path_;
};

} // namespace

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

namespace {

std::string file_digest(const fs::path &p) {
    const auto data = read_text(p);
    std::uint64_t hash = 1469598103934665603ULL;
    for (const unsigned char c : data) {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << hash;
    return out.str();
}

std::string attempt_name(int n) {
    std::ostringstream out;
    out << "attempt-" << std::setfill('0') << std::setw(3) << n;
    return out.str();
}

fs::path new_attempt(const fs::path &parent) {
    fs::create_directories(parent);
    int n = 1;
    while (fs::exists(parent / attempt_name(n))) ++n;
    const auto path = parent / attempt_name(n);
    fs::create_directories(path);
    return path;
}

void validate_pw_output(const fs::path &output) {
    const auto text = read_text(output);
    if (text.find("JOB DONE") == std::string::npos)
        throw std::runtime_error("JOB DONE not found: " + output.string());
    if (std::regex_search(text, std::regex(R"(convergence\s+NOT\s+achieved)", std::regex::icase)))
        throw std::runtime_error("SCF convergence was not achieved: " + output.string());
    if (std::regex_search(text, std::regex(R"(Error\s+in\s+routine)", std::regex::icase)))
        throw std::runtime_error("QE reported Error in routine: " + output.string());
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

std::string electronic_structure_text(const Settings &s, const fs::path &output) {
    std::ostringstream data;
    data << "electronic_energy_hartree=" << std::scientific << std::setprecision(16)
         << read_total_energy_hartree(output) << "\n"
         << "multiplicity=" << s.multiplicity << "\n"
         << "source='" << output.filename().string() << "'\n";
    return data.str();
}

fs::path density_file(const fs::path &attempt, const std::string &prefix) {
    const auto save = attempt / "out" / (prefix + ".save");
    std::vector<fs::path> found;
    for (const auto &name : {"charge-density.dat", "charge-density.hdf5"})
        if (fs::is_regular_file(save / name) && fs::file_size(save / name) > 0) found.push_back(save / name);
    if (found.size() != 1)
        throw std::runtime_error("Expected exactly one non-empty charge-density.dat or charge-density.hdf5 in " + save.string());
    return found.front();
}

std::vector<int> selected_atoms(const Settings &s, int nat) {
    auto selected = s.selected;
    if (s.system_type == "gas") {
        if (!s.multiplicity_explicit || s.multiplicity < 1)
            throw std::runtime_error("gas requires an explicit positive multiplicity in fdvib.in");
        if (!s.selected_all) throw std::runtime_error("gas requires selected_atoms='all'");
        selected.resize(nat);
        std::iota(selected.begin(), selected.end(), 1);
    }
    if (selected.empty()) throw std::runtime_error("selected_atoms is empty");
    std::set<int> unique;
    for (const int atom : selected)
        if (atom < 1 || atom > nat || !unique.insert(atom).second)
            throw std::runtime_error("Invalid/duplicate selected atom");
    return selected;
}

std::string dataset_description(const Settings &s, const QEInput &q,
                                const std::vector<int> &selected) {
    std::ostringstream out;
    out << "format=1\n"
        << "scf_digest=" << file_digest(s.scf_input) << '\n'
        << "system_type=" << s.system_type << '\n'
        << "displacement_angstrom=" << std::setprecision(17) << s.displacement << '\n'
        << "multiplicity=" << s.multiplicity << '\n'
        << "prefix=" << s.output_prefix << '\n'
        << "qe_prefix=" << q.prefix << "\nselected_atoms=";
    for (const int atom : selected) out << atom << ',';
    out << '\n';
    return out.str();
}

void initialize_dataset(const Settings &s, const QEInput &q,
                        const std::vector<int> &selected) {
    const auto state = s.workdir / "state";
    const auto wanted = dataset_description(s, q, selected);
    if (fs::exists(s.workdir) && !fs::exists(state / "dataset.state") &&
        fs::directory_iterator(s.workdir) != fs::directory_iterator())
        throw std::runtime_error("Refusing to use non-empty outdir without FDVIB state metadata: " + s.workdir.string());
    fs::create_directories(state);
    const auto dataset = state / "dataset.state";
    if (fs::exists(dataset) && read_text(dataset) != wanted)
        throw std::runtime_error("Dataset differs from the existing calculation; use a different outdir");
    if (!fs::exists(dataset)) write_text(dataset, wanted);
    const auto scf_reference = s.workdir / "scf.in.reference";
    if (!fs::exists(scf_reference)) fs::copy_file(s.scf_input, scf_reference);
    else if (read_text(scf_reference) != read_text(s.scf_input))
        throw std::runtime_error("scf.in differs from the existing calculation");
    const auto config_reference = s.workdir / "fdvib.in.reference";
    std::ostringstream cfg;
    cfg << "&FDVIB\n  scf_input='scf.in',\n  outdir='fdvib',\n"
        << "  system_type='" << s.system_type << "',\n  selected_atoms='";
    if (s.selected_all) cfg << "all";
    else for (std::size_t i = 0; i < selected.size(); ++i) cfg << (i ? "," : "") << selected[i];
    cfg << "',\n  displacement_angstrom=" << std::setprecision(17) << s.displacement << ",\n"
        << "  multiplicity=" << s.multiplicity << ",\n  prefix='" << s.output_prefix << "',\n"
        << "  pw_command='pw.x',\n  run_dynmat=.false.,\n  dynmat_command='dynmat.x',\n/\n";
    if (!fs::exists(config_reference)) write_text(config_reference, cfg.str());
    else if (read_text(config_reference) != cfg.str())
        throw std::runtime_error("FDVIB dataset snapshot is missing or modified: " + config_reference.string());
}

fs::path ensure_reference(const Settings &s, const QEInput &q) {
    const auto parent = s.workdir / "reference";
    const auto marker = parent / "complete.state";
    if (fs::exists(marker)) {
        std::istringstream in(read_text(marker));
        std::string attempt_name_saved, density_name, digest;
        in >> attempt_name_saved >> density_name >> digest;
        const auto attempt = parent / attempt_name_saved;
        const auto density = attempt / "out" / (q.prefix + ".save") / density_name;
        validate_pw_output(attempt / "scf.out");
        if (!fs::is_regular_file(density) || file_digest(density) != digest)
            throw std::runtime_error("Completed reference charge density is missing or modified: " + density.string());
        const auto electronic = s.workdir / "electronic_structure.dat";
        const auto wanted_electronic = electronic_structure_text(s, attempt / "scf.out");
        if (!fs::is_regular_file(electronic) || read_text(electronic) != wanted_electronic)
            throw std::runtime_error("Reference electronic-structure metadata is missing or modified: " + electronic.string());
        std::cout << "Preserved completed reference SCF\n";
        return density;
    }
    const auto attempt = new_attempt(parent);
    const auto input = attempt / "scf.in";
    const auto output = attempt / "scf.out";
    const auto relout = fs::relative(attempt / "out", s.root).generic_string();
    write_text(input, reference_input(q, "./" + relout));
    const auto cmd = s.pw_command + " -in " + shell_quote(fs::relative(input, s.root).string());
    std::cout << "Running unperturbed reference SCF\n";
    const int rc = shell_run(cmd, s.root, fs::relative(output, s.root));
    if (rc != 0) throw std::runtime_error("Reference SCF failed with exit code " + std::to_string(rc));
    validate_pw_output(output);
    (void)parse_forces(output, q.nat);
    const auto density = density_file(attempt, q.prefix);
    write_text(s.workdir / "electronic_structure.dat", electronic_structure_text(s, output));
    write_text(marker, attempt.filename().string() + " " + density.filename().string() + " " + file_digest(density) + "\n");
    return density;
}

void run_displacements(const Settings &s, const QEInput &q,
                       const std::vector<int> &selected, const fs::path &reference_density) {
    int completed = 0, preserved = 0;
    const auto jobs = s.workdir / "jobs";
    for (const int atom1 : selected) for (int axis = 0; axis < 3; ++axis) for (const int sign : {1, -1}) {
        const auto id = job_name(atom1, axis, sign);
        const auto parent = jobs / id;
        const auto marker = parent / "complete.state";
        const auto forces = parent / "forces.dat";
        if (fs::exists(marker)) {
            std::istringstream in(read_text(marker));
            std::string saved_attempt, digest;
            in >> saved_attempt >> digest;
            const auto output = parent / saved_attempt / "pw.out";
            const auto parsed = parse_forces(output, q.nat);
            if (!fs::is_regular_file(forces) || file_digest(forces) != digest)
                throw std::runtime_error("Completed force data is missing or modified: " + forces.string());
            (void)parsed;
            ++preserved;
            continue;
        }
        const auto attempt = new_attempt(parent);
        const auto input = attempt / "pw.in";
        const auto output = attempt / "pw.out";
        const auto relout = fs::relative(attempt / "out", s.root).generic_string();
        write_text(input, displaced_input(q, atom1 - 1, axis, sign * s.displacement, "./" + relout));
        const auto seeded = attempt / "out" / (q.prefix + ".save") / reference_density.filename();
        fs::create_directories(seeded.parent_path());
        fs::copy_file(reference_density, seeded);
        if (file_digest(seeded) != file_digest(reference_density))
            throw std::runtime_error("Copied reference charge density failed verification: " + seeded.string());
        const auto paw_src = reference_density.parent_path() / "paw.txt";
        if (fs::is_regular_file(paw_src))
            fs::copy_file(paw_src, seeded.parent_path() / "paw.txt");
        const auto cmd = s.pw_command + " -in " + shell_quote(fs::relative(input, s.root).string());
        std::cout << "Running " << id << std::endl;
        const int rc = shell_run(cmd, s.root, fs::relative(output, s.root));
        if (rc != 0) throw std::runtime_error(id + " failed with exit code " + std::to_string(rc));
        const auto parsed = parse_forces(output, q.nat);
        write_forces(forces, parsed, output);
        write_text(marker, attempt.filename().string() + " " + file_digest(forces) + "\n");
        ++completed;
    }
    std::cout << "Completed " << completed << ", preserved " << preserved << " displacement jobs\n";
}

void ensure_analysis(const Settings &s) {
    const auto marker = s.workdir / "state" / "analyze.complete";
    const auto results = s.workdir / "results";
    const auto dyn = results / (s.output_prefix + ".dynG");
    const auto result_config = results / "fdvib.in.reference";
    const auto result_electronic = results / "electronic_structure.dat";
    if (fs::exists(marker)) {
        if (!fs::is_regular_file(dyn) || !fs::is_regular_file(results / "dynmat.in"))
            throw std::runtime_error("Completed Hessian results are missing");
        std::istringstream in(read_text(marker));
        std::string dyn_digest, input_digest, config_digest, electronic_digest;
        in >> dyn_digest >> input_digest >> config_digest >> electronic_digest;
        if (file_digest(dyn) != dyn_digest || file_digest(results / "dynmat.in") != input_digest ||
            !fs::is_regular_file(result_config) || file_digest(result_config) != config_digest ||
            !fs::is_regular_file(result_electronic) || file_digest(result_electronic) != electronic_digest)
            throw std::runtime_error("Completed Hessian results were modified");
        std::cout << "Preserved completed Hessian analysis\n";
        return;
    }
    if (fs::exists(results)) {
        bool recovered = false;
        if (fs::is_regular_file(dyn) && fs::is_regular_file(results / "dynmat.in") &&
            fs::is_regular_file(result_config) && fs::is_regular_file(result_electronic)) {
            try {
                (void)read_dyn_geometry(dyn);
                const auto input = Config::load(results / "dynmat.in");
                if (input.get("fildyn") != dyn.filename().string())
                    throw std::runtime_error("dynmat.in references a different dynamical matrix");
                if (read_text(result_config) != read_text(s.workdir / "fdvib.in.reference"))
                    throw std::runtime_error("Result dataset snapshot differs");
                if (read_text(result_electronic) != read_text(s.workdir / "electronic_structure.dat"))
                    throw std::runtime_error("Result electronic-structure metadata differs");
                write_text(marker, file_digest(dyn) + " " + file_digest(results / "dynmat.in") +
                                   " " + file_digest(result_config) + " " + file_digest(result_electronic) + "\n");
                std::cout << "Recovered completed Hessian analysis\n";
                recovered = true;
            } catch (const std::exception &) {
                recovered = false;
            }
        }
        if (recovered) return;
        const auto failed_root = s.workdir / "failed";
        fs::create_directories(failed_root);
        int n = 1;
        fs::path failed;
        do { failed = failed_root / ("analysis-" + attempt_name(n++)); } while (fs::exists(failed));
        fs::rename(results, failed);
        std::cout << "Preserved incomplete Hessian results in " << failed << '\n';
    }
    analyze(s);
    write_text(marker, file_digest(dyn) + " " + file_digest(results / "dynmat.in") +
                       " " + file_digest(result_config) + " " + file_digest(result_electronic) + "\n");
}

void ensure_dynmat(const Settings &s) {
    if (!s.run_dynmat) {
        std::cout << "dynmat.x was not requested (run_dynmat=.false.)\n";
        return;
    }
    const auto state_marker = s.workdir / "state" / "dynmat.complete";
    const auto results = s.workdir / "results";
    const auto final_output = results / "dynmat.out";
    const auto final_freq = results / (s.output_prefix + ".freq.out");
    const auto dyn = results / (s.output_prefix + ".dynG");
    const auto validate_modes = [&] {
        const auto geometry = read_dyn_geometry(dyn);
        (void)parse_modes(final_freq, static_cast<int>(geometry.masses.size()));
    };
    if (fs::exists(state_marker)) {
        validate_pw_output(final_output);
        if (!fs::is_regular_file(final_freq)) throw std::runtime_error("Completed dynmat frequency output is missing");
        std::istringstream in(read_text(state_marker));
        std::string output_digest, freq_digest;
        in >> output_digest >> freq_digest;
        if (file_digest(final_output) != output_digest || file_digest(final_freq) != freq_digest)
            throw std::runtime_error("Completed dynmat results were modified");
        validate_modes();
        std::cout << "Preserved completed dynmat.x result\n";
        return;
    }
    if (fs::exists(final_output) || fs::exists(final_freq)) {
        if (fs::is_regular_file(final_output) && fs::is_regular_file(final_freq)) {
            try {
                validate_pw_output(final_output);
                validate_modes();
                write_text(state_marker, file_digest(final_output) + " " + file_digest(final_freq) + "\n");
                std::cout << "Recovered completed dynmat.x result\n";
                return;
            } catch (const std::exception &) {
                // Preserve the uncommitted files below and create a new attempt.
            }
        }
        const auto failed = new_attempt(s.workdir / "failed" / "dynmat-publish");
        if (fs::exists(final_output)) fs::rename(final_output, failed / final_output.filename());
        if (fs::exists(final_freq)) fs::rename(final_freq, failed / final_freq.filename());
        std::cout << "Preserved incomplete dynmat results in " << failed << '\n';
    }
    const auto attempt = new_attempt(s.workdir / "dynmat");
    fs::copy_file(results / (s.output_prefix + ".dynG"), attempt / (s.output_prefix + ".dynG"));
    fs::copy_file(results / "dynmat.in", attempt / "dynmat.in");
    std::cout << "Running dynmat.x\n";
    const int rc = shell_run(s.dynmat_command + " -in dynmat.in", attempt, "dynmat.out");
    if (rc != 0) throw std::runtime_error("dynmat.x failed with exit code " + std::to_string(rc));
    validate_pw_output(attempt / "dynmat.out");
    const auto freq = attempt / (s.output_prefix + ".freq.out");
    if (!fs::is_regular_file(freq) || fs::file_size(freq) == 0)
        throw std::runtime_error("dynmat.x did not produce " + freq.filename().string());
    const auto geometry = read_dyn_geometry(attempt / (s.output_prefix + ".dynG"));
    (void)parse_modes(freq, static_cast<int>(geometry.masses.size()));
    if (fs::exists(final_output) || fs::exists(final_freq))
        throw std::runtime_error("Refuse to overwrite existing dynmat result");
    fs::copy_file(attempt / "dynmat.out", final_output);
    fs::copy_file(freq, final_freq);
    write_text(state_marker, file_digest(final_output) + " " + file_digest(final_freq) + "\n");
}

} // namespace

void calculate(const Settings &s) {
    CalculationLock lock(s.workdir);
    const auto q = parse_qe_input(s.scf_input);
    const auto selected = selected_atoms(s, q.nat);
    if (s.system_type == "gas") {
        std::smatch match;
        int nspin = 1;
        if (std::regex_search(q.clean_text, match, std::regex(R"(\bnspin\s*=\s*(\d+))", std::regex::icase)))
            nspin = std::stoi(match[1]);
        double magnetization = 0.0;
        const bool has_magnetization = std::regex_search(
            q.clean_text, match, std::regex(R"(\btot_magnetization\s*=\s*([-+0-9.EeDd]+))", std::regex::icase));
        if (has_magnetization) magnetization = number(match[1]);
        if (s.multiplicity == 1) {
            if (nspin != 1 && !(nspin == 2 && has_magnetization && std::abs(magnetization) < 1e-8))
                throw std::runtime_error("Gas singlet requires nspin=1, or nspin=2 with tot_magnetization=0");
        } else if (nspin != 2 || !has_magnetization ||
                   std::abs(magnetization - (s.multiplicity - 1)) > 1e-8) {
            throw std::runtime_error("Gas multiplicity requires nspin=2 and tot_magnetization=multiplicity-1");
        }
    }
    initialize_dataset(s, q, selected);
    const auto density = ensure_reference(s, q);
    run_displacements(s, q, selected, density);
    ensure_analysis(s);
    ensure_dynmat(s);
}

} // namespace fdvib
