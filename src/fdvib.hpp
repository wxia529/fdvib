#pragma once

#include <array>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace fdvib {

namespace fs = std::filesystem;

// CODATA 2022. AMU_RY = 1 / (2 m_e[u]) for QE Rydberg atomic units.
inline constexpr double BOHR_TO_ANG = 0.529177210544;
inline constexpr double AMU_RY = 911.4442431390707;
inline constexpr double CM_TO_EV = 1.2398419843320026e-4;
inline constexpr double KB_EV = 8.617333262145179e-5;
inline constexpr double KB_SI = 1.380649e-23;
inline constexpr double H_SI = 6.62607015e-34;
inline constexpr double AMU_KG = 1.66053906892e-27;
inline constexpr double ATM_PA = 101325.0;
inline constexpr double PI = 3.141592653589793238462643383279502884;

using Vec3 = std::array<double, 3>;

std::string trim(const std::string &s);
std::string lower(std::string s);
std::string unquote(std::string s);
std::string strip_comment(const std::string &line);
double number(std::string s);
std::string read_text(const fs::path &p);
void write_text(const fs::path &p, const std::string &s);
bool is_standard_element(const std::string &symbol);
std::string standard_element_symbol(const std::string &label, const std::string &context);
int atomic_number_from_label(const std::string &label, const std::string &context);

struct Config {
    std::map<std::string, std::string> v;

    static Config load(const fs::path &p);
    static Config load_qe_namelist(const fs::path &p);
    bool has(const std::string &k) const;
    std::string get(const std::string &k, const std::string &d = {}) const;
    double real(const std::string &k, double d) const;
    int integer(const std::string &k, int d) const;
    void require_only(const std::set<std::string> &allowed,
                      const std::string &context) const;
};

std::vector<int> integer_list(std::string s);

struct Species { std::string symbol, pseudo; double mass{}; };
struct Atom {
    std::string symbol;
    Vec3 r{};
    Vec3 input_r{};
    std::vector<std::string> extra;
    int type{};
};

struct QEInput {
    std::string text, clean_text;
    std::string prefix{"pwscf"};
    std::vector<std::string> lines;
    int nat{}, ntyp{}, pos_header{-1}, pos_start{-1}, cell_header{-1};
    std::string positions_unit{"angstrom"}, cell_unit{"angstrom"};
    double alat_angstrom{};
    std::vector<Species> species;
    std::vector<Atom> atoms;
    std::array<Vec3, 3> cell{};
};

QEInput parse_qe_input(const fs::path &p);
std::string format_position(const Atom &a, const Vec3 &coord);
std::string reference_input(const QEInput &q, const std::string &outdir,
                            const fs::path &source_dir, const fs::path &run_dir);
std::string displaced_input(const QEInput &q, int atom, int axis, double shift,
                            const std::string &outdir, const fs::path &source_dir,
                            const fs::path &run_dir);
std::string shell_quote(const std::string &s);
int shell_run(const std::string &cmd, const fs::path &cwd, const fs::path &stdout_path);
std::string file_digest(const fs::path &p);
std::string numbered_name(const std::string &task, int n);
bool is_numbered_name(const std::string &name, const std::string &task);
fs::path new_numbered_directory(const fs::path &parent, const std::string &task);
std::vector<fs::path> numbered_directories(const fs::path &parent, const std::string &task);
void validate_qe_output(const fs::path &output);
std::vector<Vec3> parse_forces(const fs::path &output, int nat);
void write_forces(const fs::path &p, const std::vector<Vec3> &f, const fs::path &source);
std::vector<Vec3> read_forces(const fs::path &p, int nat);
double read_total_energy_hartree(const fs::path &output);

struct Settings {
    fs::path config_path, root, scf_input, workdir;
    std::string system_type, pw_command, dynmat_command, output_prefix;
    std::vector<int> selected;
    double displacement{};
    int multiplicity{};
    bool multiplicity_explicit{}, selected_all{}, run_dynmat{};
};

Settings settings(const fs::path &config_path, const fs::path &root_override = {});
std::string job_name(int atom1, int axis, int sign);

struct ResultMetadata {
    std::string program{"qe"};
    double electronic_energy_hartree{};
    int multiplicity{1};
    std::string mode_selection{"all"};
    std::string selected_atoms{"all"};
};

ResultMetadata result_metadata(const fs::path &results, bool required);

struct ResultFiles {
    fs::path dyn;
    fs::path freq;
};

ResultFiles result_files(const fs::path &results, const std::string &context);

struct Mode { double freq{}; std::vector<Vec3> displacement; };
struct DynGeometry {
    std::vector<double> masses;
    std::vector<Vec3> r_bohr;
    std::vector<std::string> symbols;
};

struct ModeSelection {
    std::vector<const Mode *> modes;
    std::string classification;
    double largest_removed{};
};

ModeSelection select_gas_internal_modes(const DynGeometry &geometry,
                                        const std::vector<Mode> &modes,
                                        const std::string &context);
ModeSelection select_shm_modes(const ResultMetadata &metadata,
                               const DynGeometry &geometry,
                               const std::vector<Mode> &modes);
ModeSelection select_fakeg_modes(const ResultMetadata &metadata,
                                 const DynGeometry &geometry,
                                 const std::vector<Mode> &modes);

std::vector<Mode> parse_modes(const fs::path &p, int nat);
DynGeometry read_dyn_geometry(const fs::path &p);

void analyze(const Settings &s);
void calculate(const Settings &s);
void modes(const fs::path &results);
void thermo(const fs::path &results, const fs::path &thermo_input);
void shm(const fs::path &results);
void fakeg(const fs::path &results);

} // namespace fdvib
