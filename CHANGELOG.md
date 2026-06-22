# Changelog

All notable changes to FDVIB are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.1] - 2026-06-22

### Changed

- Build Linux release packages with the manylinux2014 toolchain for glibc 2.17
  compatibility.
- Link `libstdc++` and `libgcc` statically in Linux release packages.
- Include the glibc compatibility baseline in Linux archive names.

## [0.1.0] - 2026-06-22

### Added

- Cartesian finite-difference displacement generation and Quantum ESPRESSO
  `pw.x` execution.
- Frozen-environment local harmonic calculations for selected atoms.
- Gas-phase rigid-rotor harmonic-oscillator thermochemistry.
- QE-compatible dynamical-matrix and `dynmat.x` input generation.
- Compact Molden normal-mode output with signed frequencies.
- Harmonic and frequency-floor treatments for positive low-frequency modes.
- Configuration snapshots and overwrite protection for calculation data.
- Reference documentation and local and gas-phase configuration examples.
- Tag-triggered Linux release packaging through GitHub Actions.
- BSD 3-Clause licensing.

[Unreleased]: https://github.com/wxia529/fdvib/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/wxia529/fdvib/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/wxia529/fdvib/releases/tag/v0.1.0
