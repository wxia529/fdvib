// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Wanting Xia

#include "fdvib.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>

namespace fdvib {

namespace {

std::array<double, 3> inertia_eigenvalues(const DynGeometry &g, const std::string &context) {
    const double total_mass = std::accumulate(g.masses.begin(), g.masses.end(), 0.0);
    if (!(total_mass > 0.0)) throw std::runtime_error("Invalid total mass for " + context);
    Vec3 center{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < g.masses.size(); ++i)
        for (int k = 0; k < 3; ++k) center[k] += g.masses[i] * g.r_bohr[i][k] / total_mass;
    Vec3 diagonal{0.0, 0.0, 0.0}, off{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < g.masses.size(); ++i) {
        const double x = g.r_bohr[i][0] - center[0];
        const double y = g.r_bohr[i][1] - center[1];
        const double z = g.r_bohr[i][2] - center[2];
        const double mass = g.masses[i];
        diagonal[0] += mass * (y*y + z*z);
        diagonal[1] += mass * (x*x + z*z);
        diagonal[2] += mass * (x*x + y*y);
        off[0] -= mass*x*y; off[1] -= mass*x*z; off[2] -= mass*y*z;
    }
    const double p1 = off[0]*off[0] + off[1]*off[1] + off[2]*off[2];
    if (p1 == 0.0) {
        std::array<double, 3> values{diagonal[0], diagonal[1], diagonal[2]};
        std::sort(values.begin(), values.end());
        return values;
    }
    const double q = (diagonal[0] + diagonal[1] + diagonal[2]) / 3.0;
    const double p2 = (diagonal[0]-q)*(diagonal[0]-q) + (diagonal[1]-q)*(diagonal[1]-q) +
                      (diagonal[2]-q)*(diagonal[2]-q) + 2.0*p1;
    const double p = std::sqrt(p2 / 6.0);
    const double a00=(diagonal[0]-q)/p, a11=(diagonal[1]-q)/p, a22=(diagonal[2]-q)/p;
    const double a01=off[0]/p, a02=off[1]/p, a12=off[2]/p;
    const double determinant = a00*a11*a22 + 2*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
    const double phi = std::acos(std::clamp(determinant/2.0, -1.0, 1.0)) / 3.0;
    const double e3=q+2*p*std::cos(phi), e1=q+2*p*std::cos(phi+2*PI/3), e2=3*q-e1-e3;
    std::array<double, 3> values{e1,e2,e3};
    std::sort(values.begin(), values.end());
    return values;
}

ModeSelection retain_largest_modes(const std::vector<Mode> &modes, std::size_t keep,
                                   const std::string &classification,
                                   const std::string &context) {
    if (keep > modes.size())
        throw std::runtime_error("Requested " + context + " vibration count exceeds available modes");
    std::vector<std::size_t> order(modes.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return std::abs(modes[a].freq) < std::abs(modes[b].freq);
    });
    std::vector<bool> removed(modes.size(), false);
    ModeSelection selected;
    selected.classification = classification;
    selected.modes.reserve(keep);
    for (std::size_t i = 0; i < modes.size() - keep; ++i) {
        removed[order[i]] = true;
        selected.largest_removed = std::max(selected.largest_removed, std::abs(modes[order[i]].freq));
    }
    for (std::size_t i = 0; i < modes.size(); ++i)
        if (!removed[i]) selected.modes.push_back(&modes[i]);
    return selected;
}

ModeSelection select_nonzero_modes(const std::vector<Mode> &modes, const std::string &classification) {
    ModeSelection selected;
    selected.classification = classification;
    selected.modes.reserve(modes.size());
    for (const auto &mode : modes)
        if (mode.freq != 0.0) selected.modes.push_back(&mode);
    return selected;
}

std::vector<int> metadata_selected_atoms(const ResultMetadata &metadata, std::size_t nat) {
    std::vector<int> selected;
    if (metadata.selected_atoms == "all") {
        selected.resize(nat);
        std::iota(selected.begin(), selected.end(), 1);
    } else {
        selected = integer_list(metadata.selected_atoms);
    }
    if (selected.empty()) throw std::runtime_error("Local SHM export requires selected_atoms in metadata.dat");
    std::set<int> unique;
    for (const int atom : selected)
        if (atom < 1 || static_cast<std::size_t>(atom) > nat || !unique.insert(atom).second)
            throw std::runtime_error("Invalid/duplicate selected atom in metadata.dat");
    return selected;
}

} // namespace

ModeSelection select_gas_internal_modes(const DynGeometry &geometry,
                                        const std::vector<Mode> &modes,
                                        const std::string &context) {
    const auto nat = geometry.masses.size();
    if (nat == 1) {
        ModeSelection selected;
        selected.classification = "atom";
        return selected;
    }
    const auto inertia = inertia_eigenvalues(geometry, context);
    const bool linear = std::any_of(inertia.begin(), inertia.end(), [](double value) { return value < 0.001; });
    return retain_largest_modes(modes, 3*nat - (linear ? 5 : 6),
                                linear ? "linear" : "nonlinear", context);
}

ModeSelection select_shm_modes(const ResultMetadata &metadata,
                               const DynGeometry &geometry,
                               const std::vector<Mode> &modes) {
    if (metadata.mode_selection == "gas")
        return select_gas_internal_modes(geometry, modes, "SHM export");
    if (metadata.mode_selection == "local")
        return retain_largest_modes(modes, 3 * metadata_selected_atoms(metadata, geometry.masses.size()).size(),
                                    "local", "SHM export");
    if (metadata.mode_selection == "all")
        return select_nonzero_modes(modes, "all");
    throw std::runtime_error("mode_selection must be all, gas, or local in metadata.dat");
}

ModeSelection select_fakeg_modes(const ResultMetadata &metadata,
                                 const DynGeometry &geometry,
                                 const std::vector<Mode> &modes) {
    if (metadata.mode_selection == "gas")
        return select_gas_internal_modes(geometry, modes, "fake Gaussian export");
    if (metadata.mode_selection == "local" || metadata.mode_selection == "all")
        return select_nonzero_modes(modes, metadata.mode_selection);
    throw std::runtime_error("mode_selection must be all, gas, or local in metadata.dat");
}

} // namespace fdvib
