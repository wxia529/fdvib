#include "fdvib.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace fdvib {

std::string file_digest(const fs::path &p) {
    std::ifstream input(p, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read " + p.string());
    std::uint64_t hash = 1469598103934665603ULL;
    std::array<char, 65536> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = input.gcount();
        for (std::streamsize i = 0; i < count; ++i) {
            hash ^= static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]);
            hash *= 1099511628211ULL;
        }
    }
    if (!input.eof()) throw std::runtime_error("Failed while reading " + p.string());
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << hash;
    return out.str();
}

std::string numbered_name(const std::string &task, int n) {
    std::ostringstream out;
    out << task << '_' << std::setfill('0') << std::setw(3) << n;
    return out.str();
}

bool is_numbered_name(const std::string &name, const std::string &task) {
    const auto prefix = task + '_';
    if (name.rfind(prefix, 0) != 0 || name.size() < prefix.size() + 3) return false;
    return std::all_of(name.begin() + static_cast<std::ptrdiff_t>(prefix.size()), name.end(),
                       [](unsigned char c) { return std::isdigit(c) != 0; });
}

fs::path new_numbered_directory(const fs::path &parent, const std::string &task) {
    fs::create_directories(parent);
    int n = 1;
    while (fs::exists(parent / numbered_name(task, n))) ++n;
    const auto path = parent / numbered_name(task, n);
    fs::create_directories(path);
    return path;
}

std::vector<fs::path> numbered_directories(const fs::path &parent, const std::string &task) {
    std::vector<fs::path> found;
    if (!fs::is_directory(parent)) return found;
    for (const auto &entry : fs::directory_iterator(parent))
        if (entry.is_directory() && is_numbered_name(entry.path().filename().string(), task))
            found.push_back(entry.path());
    std::sort(found.rbegin(), found.rend());
    return found;
}

} // namespace fdvib
