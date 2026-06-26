#include "fdvib.hpp"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
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

namespace {

std::string metadata_text(const Settings &s, const fs::path &output) {
    std::ostringstream data;
    data << "program = qe\n"
         << "electronic_energy_hartree = " << std::scientific << std::setprecision(16)
         << read_total_energy_hartree(output) << "\n"
         << "multiplicity = " << s.multiplicity << "\n"
         << "mode_selection = " << s.system_type << "\n"
         << "selected_atoms = ";
    if (s.selected_all) data << "all";
    else for (std::size_t i = 0; i < s.selected.size(); ++i) data << (i ? "," : "") << s.selected[i];
    data << "\n";
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
    out << "format=2\n"
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
    cfg << "scf_input = scf.in\noutdir = fdvib\n"
        << "system_type = " << s.system_type << "\nselected_atoms = ";
    if (s.selected_all) cfg << "all";
    else for (std::size_t i = 0; i < selected.size(); ++i) cfg << (i ? "," : "") << selected[i];
    cfg << "\ndisplacement_angstrom = " << std::setprecision(17) << s.displacement << "\n"
        << "multiplicity = " << s.multiplicity << "\nprefix = " << s.output_prefix << "\n"
        << "pw_command = pw.x\nrun_dynmat = false\ndynmat_command = dynmat.x\n";
    if (!fs::exists(config_reference)) write_text(config_reference, cfg.str());
    else if (read_text(config_reference) != cfg.str())
        throw std::runtime_error("FDVIB dataset snapshot is missing or modified: " + config_reference.string());
}

struct ReferenceSeed {
    fs::path density;
    fs::path paw;
    std::string density_digest;
    std::string paw_digest;
};

ReferenceSeed ensure_reference(const Settings &s, const QEInput &q) {
    const auto calculations = s.workdir / "calculations";
    const auto marker = s.workdir / "state" / "init_scf.complete";
    if (fs::exists(marker)) {
        std::istringstream in(read_text(marker));
        std::string attempt_name_saved, density_name, density_digest, paw_name, paw_digest, extra;
        if (!(in >> attempt_name_saved >> density_name >> density_digest >> paw_name >> paw_digest) || in >> extra ||
            !is_numbered_name(attempt_name_saved, "init_scf"))
            throw std::runtime_error("Invalid reference completion snapshot: " + marker.string());
        const auto attempt = calculations / attempt_name_saved;
        const auto density = attempt / "out" / (q.prefix + ".save") / density_name;
        const auto paw = density.parent_path() / "paw.txt";
        validate_qe_output(attempt / "scf.out");
        if (!fs::is_regular_file(density) || file_digest(density) != density_digest)
            throw std::runtime_error("Completed reference charge density is missing or modified: " + density.string());
        if (paw_name == "-") {
            if (paw_digest != "-" || fs::exists(paw))
                throw std::runtime_error("Completed reference PAW snapshot is inconsistent: " + paw.string());
        } else if (paw_name != "paw.txt" || !fs::is_regular_file(paw) || fs::file_size(paw) == 0 ||
                   file_digest(paw) != paw_digest) {
            throw std::runtime_error("Completed reference PAW data is missing or modified: " + paw.string());
        }
        const auto metadata = s.workdir / "metadata.dat";
        const auto wanted_metadata = metadata_text(s, attempt / "scf.out");
        if (!fs::is_regular_file(metadata) || read_text(metadata) != wanted_metadata)
            throw std::runtime_error("Reference metadata is missing or modified: " + metadata.string());
        std::cout << "Preserved completed reference SCF\n";
        return {density, paw_name == "paw.txt" ? paw : fs::path{}, density_digest,
                paw_name == "paw.txt" ? paw_digest : std::string{}};
    }
    const auto commit = [&](const fs::path &calculation, bool recovered) {
        const auto output = calculation / "scf.out";
        validate_qe_output(output);
        (void)parse_forces(output, q.nat);
        const auto density = density_file(calculation, q.prefix);
        const auto paw = density.parent_path() / "paw.txt";
        if (fs::exists(paw) && (!fs::is_regular_file(paw) || fs::file_size(paw) == 0))
            throw std::runtime_error("Reference PAW data is not a non-empty regular file: " + paw.string());
        write_text(s.workdir / "metadata.dat", metadata_text(s, output));
        const bool has_paw = fs::is_regular_file(paw);
        const auto density_digest = file_digest(density);
        const auto paw_digest = has_paw ? file_digest(paw) : std::string{};
        write_text(marker, calculation.filename().string() + " " + density.filename().string() + " " +
                           density_digest + " " + (has_paw ? "paw.txt " + paw_digest : "- -") + "\n");
        if (recovered) std::cout << "Recovered completed initial SCF from " << calculation.filename().string() << '\n';
        return ReferenceSeed{density, has_paw ? paw : fs::path{}, density_digest, paw_digest};
    };
    for (const auto &calculation : numbered_directories(calculations, "init_scf")) {
        try {
            return commit(calculation, true);
        } catch (const std::exception &) {
            // Retain incomplete calculations and try a fresh numbered directory below.
        }
    }
    const auto attempt = new_numbered_directory(calculations, "init_scf");
    const auto input = attempt / "scf.in";
    write_text(input, reference_input(q, "./out", s.root, attempt));
    const auto cmd = s.pw_command + " -inp scf.in";
    std::cout << "Running unperturbed reference SCF\n";
    const int rc = shell_run(cmd, attempt, "scf.out");
    if (rc != 0) throw std::runtime_error("Reference SCF failed with exit code " + std::to_string(rc));
    return commit(attempt, false);
}

void run_displacements(const Settings &s, const QEInput &q,
                       const std::vector<int> &selected, const ReferenceSeed &reference) {
    int completed = 0, preserved = 0;
    const auto calculations = s.workdir / "calculations";
    for (const int atom1 : selected) for (int axis = 0; axis < 3; ++axis) for (const int sign : {1, -1}) {
        const auto id = job_name(atom1, axis, sign);
        const auto marker = s.workdir / "state" / (id + ".complete");
        if (fs::exists(marker)) {
            std::istringstream in(read_text(marker));
            std::string saved_attempt, digest, extra;
            if (!(in >> saved_attempt >> digest) || in >> extra || !is_numbered_name(saved_attempt, id))
                throw std::runtime_error("Invalid displacement completion snapshot: " + marker.string());
            const auto attempt = calculations / saved_attempt;
            const auto output = attempt / "pw.out";
            const auto forces = attempt / "forces.dat";
            const auto parsed = parse_forces(output, q.nat);
            if (!fs::is_regular_file(forces) || file_digest(forces) != digest)
                throw std::runtime_error("Completed force data is missing or modified: " + forces.string());
            (void)parsed;
            ++preserved;
            continue;
        }
        const auto commit = [&](const fs::path &calculation, bool recovered) {
            const auto output = calculation / "pw.out";
            const auto forces = calculation / "forces.dat";
            const auto parsed = parse_forces(output, q.nat);
            write_forces(forces, parsed, output);
            write_text(marker, calculation.filename().string() + " " + file_digest(forces) + "\n");
            if (recovered) std::cout << "Recovered completed " << id << " from " << calculation.filename().string() << '\n';
        };
        bool recovered = false;
        for (const auto &calculation : numbered_directories(calculations, id)) {
            try {
                commit(calculation, true);
                recovered = true;
                break;
            } catch (const std::exception &) {
                // Retain incomplete calculations and try a fresh numbered directory below.
            }
        }
        if (recovered) {
            ++completed;
            continue;
        }
        const auto attempt = new_numbered_directory(calculations, id);
        const auto input = attempt / "pw.in";
        write_text(input, displaced_input(q, atom1 - 1, axis, sign * s.displacement,
                                          "./out", s.root, attempt));
        const auto seeded = attempt / "out" / (q.prefix + ".save") / reference.density.filename();
        fs::create_directories(seeded.parent_path());
        fs::copy_file(reference.density, seeded);
        if (file_digest(seeded) != reference.density_digest)
            throw std::runtime_error("Copied reference charge density failed verification: " + seeded.string());
        if (!reference.paw.empty()) {
            const auto seeded_paw = seeded.parent_path() / "paw.txt";
            fs::copy_file(reference.paw, seeded_paw);
            if (file_digest(seeded_paw) != reference.paw_digest)
                throw std::runtime_error("Copied reference PAW data failed verification: " + seeded_paw.string());
        }
        const auto cmd = s.pw_command + " -inp pw.in";
        std::cout << "Running " << id << std::endl;
        const int rc = shell_run(cmd, attempt, "pw.out");
        if (rc != 0) throw std::runtime_error(id + " failed with exit code " + std::to_string(rc));
        commit(attempt, false);
        ++completed;
    }
    std::cout << "Completed " << completed << ", preserved " << preserved << " displacement jobs\n";
}

void ensure_analysis(const Settings &s) {
    const auto marker = s.workdir / "state" / "analyze.complete";
    const auto results = s.workdir / "results";
    const auto dyn = results / (s.output_prefix + ".dynG");
    const auto result_metadata = results / "metadata.dat";
    if (fs::exists(marker)) {
        if (!fs::is_regular_file(dyn) || !fs::is_regular_file(results / "dynmat.in"))
            throw std::runtime_error("Completed Hessian results are missing");
        std::istringstream in(read_text(marker));
        std::string dyn_digest, input_digest, metadata_digest, extra;
        if (!(in >> dyn_digest >> input_digest >> metadata_digest) || in >> extra)
            throw std::runtime_error("Invalid Hessian completion snapshot: " + marker.string());
        if (file_digest(dyn) != dyn_digest || file_digest(results / "dynmat.in") != input_digest ||
            !fs::is_regular_file(result_metadata) || file_digest(result_metadata) != metadata_digest)
            throw std::runtime_error("Completed Hessian results were modified");
        std::cout << "Preserved completed Hessian analysis\n";
        return;
    }
    if (fs::exists(results)) {
        bool recovered = false;
        if (fs::is_regular_file(dyn) && fs::is_regular_file(results / "dynmat.in") &&
            fs::is_regular_file(result_metadata)) {
            try {
                (void)read_dyn_geometry(dyn);
                const auto input = Config::load_qe_namelist(results / "dynmat.in");
                if (input.get("fildyn") != dyn.filename().string())
                    throw std::runtime_error("dynmat.in references a different dynamical matrix");
                if (read_text(result_metadata) != read_text(s.workdir / "metadata.dat"))
                    throw std::runtime_error("Result metadata differs");
                write_text(marker, file_digest(dyn) + " " + file_digest(results / "dynmat.in") +
                                   " " + file_digest(result_metadata) + "\n");
                std::cout << "Recovered completed Hessian analysis\n";
                recovered = true;
            } catch (const std::exception &) {
                recovered = false;
            }
        }
        if (recovered) return;
        const auto failed_root = s.workdir / "failed";
        const auto failed = new_numbered_directory(failed_root, "analysis");
        fs::rename(results, failed);
        std::cout << "Preserved incomplete Hessian results in " << failed << '\n';
    }
    analyze(s);
    write_text(marker, file_digest(dyn) + " " + file_digest(results / "dynmat.in") +
                       " " + file_digest(result_metadata) + "\n");
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
    const auto calculations = s.workdir / "calculations";
    const auto validate_modes = [&] {
        const auto geometry = read_dyn_geometry(dyn);
        (void)parse_modes(final_freq, static_cast<int>(geometry.masses.size()));
    };
    const auto validate_calculation = [&](const fs::path &calculation) {
        const auto calculation_dyn = calculation / dyn.filename();
        const auto calculation_input = calculation / "dynmat.in";
        const auto calculation_output = calculation / "dynmat.out";
        const auto calculation_freq = calculation / final_freq.filename();
        if (!fs::is_regular_file(calculation_dyn) || !fs::is_regular_file(calculation_input) ||
            file_digest(calculation_dyn) != file_digest(dyn) ||
            file_digest(calculation_input) != file_digest(results / "dynmat.in"))
            throw std::runtime_error("dynmat calculation inputs differ from Hessian results: " + calculation.string());
        validate_qe_output(calculation_output);
        if (!fs::is_regular_file(calculation_freq) || fs::file_size(calculation_freq) == 0)
            throw std::runtime_error("dynmat calculation frequency output is missing: " + calculation_freq.string());
        const auto geometry = read_dyn_geometry(calculation_dyn);
        (void)parse_modes(calculation_freq, static_cast<int>(geometry.masses.size()));
    };
    if (fs::exists(state_marker)) {
        std::istringstream in(read_text(state_marker));
        std::string calculation_name, output_digest, freq_digest, extra;
        if (!(in >> calculation_name >> output_digest >> freq_digest) || in >> extra ||
            !is_numbered_name(calculation_name, "dynmat"))
            throw std::runtime_error("Invalid dynmat completion snapshot: " + state_marker.string());
        const auto calculation = calculations / calculation_name;
        validate_calculation(calculation);
        validate_qe_output(final_output);
        if (!fs::is_regular_file(final_freq)) throw std::runtime_error("Completed dynmat frequency output is missing");
        if (file_digest(final_output) != output_digest || file_digest(final_freq) != freq_digest)
            throw std::runtime_error("Completed dynmat results were modified");
        if (file_digest(calculation / "dynmat.out") != output_digest ||
            file_digest(calculation / final_freq.filename()) != freq_digest)
            throw std::runtime_error("Completed dynmat calculation was modified: " + calculation.string());
        validate_modes();
        std::cout << "Preserved completed dynmat.x result\n";
        return;
    }
    if (fs::exists(final_output) || fs::exists(final_freq)) {
        if (fs::is_regular_file(final_output) && fs::is_regular_file(final_freq)) {
            for (const auto &calculation : numbered_directories(calculations, "dynmat")) {
                try {
                    validate_qe_output(final_output);
                    validate_modes();
                    validate_calculation(calculation);
                    if (file_digest(final_output) != file_digest(calculation / "dynmat.out") ||
                        file_digest(final_freq) != file_digest(calculation / final_freq.filename()))
                        continue;
                    write_text(state_marker, calculation.filename().string() + " " +
                               file_digest(final_output) + " " + file_digest(final_freq) + "\n");
                    std::cout << "Recovered completed dynmat.x result from " << calculation.filename().string() << '\n';
                    return;
                } catch (const std::exception &) {
                    // Try another retained dynmat calculation.
                }
            }
        }
        const auto failed = new_numbered_directory(s.workdir / "failed", "dynmat_publish");
        if (fs::exists(final_output)) fs::rename(final_output, failed / final_output.filename());
        if (fs::exists(final_freq)) fs::rename(final_freq, failed / final_freq.filename());
        std::cout << "Preserved incomplete dynmat results in " << failed << '\n';
    }
    for (const auto &calculation : numbered_directories(calculations, "dynmat")) {
        try {
            validate_calculation(calculation);
            fs::copy_file(calculation / "dynmat.out", final_output);
            fs::copy_file(calculation / final_freq.filename(), final_freq);
            write_text(state_marker, calculation.filename().string() + " " +
                       file_digest(final_output) + " " + file_digest(final_freq) + "\n");
            std::cout << "Recovered completed dynmat.x calculation from " << calculation.filename().string() << '\n';
            return;
        } catch (const std::exception &) {
            // Retain incomplete calculations and create a fresh directory below.
        }
    }
    const auto attempt = new_numbered_directory(calculations, "dynmat");
    fs::copy_file(results / (s.output_prefix + ".dynG"), attempt / (s.output_prefix + ".dynG"));
    fs::copy_file(results / "dynmat.in", attempt / "dynmat.in");
    std::cout << "Running dynmat.x\n";
    const int rc = shell_run(s.dynmat_command + " -inp dynmat.in", attempt, "dynmat.out");
    if (rc != 0) throw std::runtime_error("dynmat.x failed with exit code " + std::to_string(rc));
    validate_qe_output(attempt / "dynmat.out");
    const auto freq = attempt / (s.output_prefix + ".freq.out");
    if (!fs::is_regular_file(freq) || fs::file_size(freq) == 0)
        throw std::runtime_error("dynmat.x did not produce " + freq.filename().string());
    const auto geometry = read_dyn_geometry(attempt / (s.output_prefix + ".dynG"));
    (void)parse_modes(freq, static_cast<int>(geometry.masses.size()));
    if (fs::exists(final_output) || fs::exists(final_freq))
        throw std::runtime_error("Refuse to overwrite existing dynmat result");
    fs::copy_file(attempt / "dynmat.out", final_output);
    fs::copy_file(freq, final_freq);
    write_text(state_marker, attempt.filename().string() + " " +
               file_digest(final_output) + " " + file_digest(final_freq) + "\n");
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
    const auto reference = ensure_reference(s, q);
    run_displacements(s, q, selected, reference);
    ensure_analysis(s);
    ensure_dynmat(s);
}

} // namespace fdvib
