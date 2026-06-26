#include "fdvib.hpp"

#include <cmath>
#include <stdexcept>

namespace fdvib {

namespace {

std::string count_message(const std::vector<fs::path> &paths) {
    if (paths.empty()) return "0";
    std::string out = std::to_string(paths.size()) + " (";
    for (std::size_t i = 0; i < paths.size(); ++i) {
        if (i) out += ", ";
        out += paths[i].filename().string();
    }
    out += ")";
    return out;
}

} // namespace

ResultMetadata result_metadata(const fs::path &results, bool required) {
    ResultMetadata metadata;
    const auto path = results / "metadata.dat";
    if (!fs::exists(path)) {
        if (required) throw std::runtime_error("metadata.dat is required in " + results.string());
        return metadata;
    }
    const auto c = Config::load(path);
    c.require_only({"program", "electronic_energy_hartree", "multiplicity",
                    "mode_selection", "selected_atoms"}, "metadata.dat");
    metadata.program = lower(c.get("program", metadata.program));
    if (metadata.program != "qe")
        throw std::runtime_error("Unsupported program in metadata.dat: " + metadata.program);
    metadata.electronic_energy_hartree =
        c.real("electronic_energy_hartree", metadata.electronic_energy_hartree);
    metadata.multiplicity = c.integer("multiplicity", metadata.multiplicity);
    metadata.mode_selection = lower(c.get("mode_selection", metadata.mode_selection));
    metadata.selected_atoms = lower(c.get("selected_atoms", metadata.selected_atoms));
    if (metadata.mode_selection != "all" && metadata.mode_selection != "gas" &&
        metadata.mode_selection != "local")
        throw std::runtime_error("mode_selection must be all, gas, or local in metadata.dat");
    if (!std::isfinite(metadata.electronic_energy_hartree) || metadata.multiplicity <= 0)
        throw std::runtime_error("Invalid metadata.dat");
    return metadata;
}

ResultFiles result_files(const fs::path &results, const std::string &context) {
    if (!fs::is_directory(results)) throw std::runtime_error("Not a result directory: " + results.string());
    std::vector<fs::path> dyns, freqs;
    for (const auto &entry : fs::directory_iterator(results)) {
        if (!entry.is_regular_file()) continue;
        const auto name = entry.path().filename().string();
        if (name.size() >= 5 && name.substr(name.size() - 5) == ".dynG") dyns.push_back(entry.path());
        if (name.size() >= 9 && name.substr(name.size() - 9) == ".freq.out") freqs.push_back(entry.path());
    }
    if (dyns.size() != 1 || freqs.size() != 1)
        throw std::runtime_error(context + " requires exactly one .dynG and one .freq.out in " +
                                 results.string() + "; found .dynG=" + count_message(dyns) +
                                 ", .freq.out=" + count_message(freqs));
    return {dyns.front(), freqs.front()};
}

} // namespace fdvib
