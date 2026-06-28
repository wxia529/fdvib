# Changelog

All notable changes to FDVIB are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- Make `fdvib modes` write a CP2K-style `<prefix>.mol`
  containing atoms, signed frequencies, normal coordinates, and periodic cell
  vectors for every mode-selection type.
- Omit the Molden `[INT]` section because FDVIB does not calculate IR
  intensities.
- Separate QE `.dynG` and `dynmat.x` output parsing from Hessian construction,
  make the QE-specific reader interfaces explicit, and represent selected
  modes by stable indices instead of non-owning pointers.

### Removed

- Remove the built-in `fdvib fakeg` Gaussian-like export in favor of the
  generated CP2K-style Molden file.
- Remove the bundled regression test suite and CTest integration.

## [0.4.3] - 2026-06-26

### Changed

- Share gas internal-mode selection between SHM and fake Gaussian exports.
- Allow post-processing exports (`modes`, `thermo`, `shm`, and `fakeg`) to be
  re-run by overwriting their generated output files.
- Split QE input handling, QE output parsing, process execution, and run-state
  helpers out of the calculation workflow implementation.

## [0.4.2] - 2026-06-26

### Changed

- Make fake Gaussian exports use gas molecular degree-of-freedom mode
  selection when `metadata.dat` has `mode_selection=gas`; `all` and `local`
  exports still omit exact zero frequencies.

## [0.4.1] - 2026-06-26

### Added

- Add `fdvib fakeg RESULTS_DIR` to export a minimal Gaussian-like
  `<prefix>_fake.out` from FDVIB `.dynG` and `.freq.out` results for
  GaussView vibration visualization.

## [0.4.0] - 2026-06-26

### Added

- Add `results/metadata.dat` as compact result metadata for SHM export,
  including `program`, `mode_selection`, electronic energy, multiplicity, and
  selected atoms.
- Add a forward-compatible SHM `program` metadata field. The interface defaults
  to `qe`; other programs are rejected until readers are implemented.

### Changed

- Replace FDVIB and thermo namelist-style input files with strict key-value
  files.
- Make `fdvib shm` depend only on one `.dynG` file and one `.freq.out` file,
  with optional `metadata.dat`.
- Default SHM mode selection to `all`, writing all nonzero frequencies and
  omitting only frequencies exactly equal to `0.0`.
- Stop using `fdvib.in.reference`, `dynmat.in`, and `electronic_structure.dat`
  as SHM export inputs.
- Stop publishing `fdvib.in.reference` and `electronic_structure.dat` in
  `results/`; `metadata.dat` is the result-level metadata interface.
- Run the release-package test suite inside the manylinux build workflow
  before creating a GitHub release.

### Fixed

- Harden QE input and dynG parsing against truncated files, non-positive atom
  counts, invalid type indices, and non-finite coordinates.
- Report the concrete `.dynG` and `.freq.out` files found when result
  directories do not contain exactly one of each.

## [0.3.4] - 2026-06-25

### Added

- Support explicit QE position units `angstrom`, `bohr`, `alat`, and `crystal`
  when generating finite-displacement inputs.
- Support explicit QE cell units `angstrom`, `bohr`, and `alat`.
- Accept QE card-unit spellings with bare units, parentheses, or braces, such
  as `ATOMIC_POSITIONS (bohr)` and `ATOMIC_POSITIONS {crystal}`.

### Changed

- Preserve the user's original position unit when writing displaced QE inputs,
  converting the Cartesian FDVIB displacement into `bohr`, `alat`, or
  fractional `crystal` increments as needed.
- Reject omitted/deprecated QE card units and `ATOMIC_POSITIONS crystal_sg`
  explicitly instead of silently assuming an unsafe interpretation.

## [0.3.3] - 2026-06-24

### Fixed

- Expand home-relative QE `pseudo_dir` values before relocating generated
  inputs, preventing paths such as `../../../~/pseudo`.

## [0.3.2] - 2026-06-24

### Fixed

- Include optional reference `paw.txt` data in completion snapshots and verify
  each copy used to seed a displaced PAW calculation.
- Run QE inside each numbered calculation directory with `outdir='./out'`,
  allowing a generated `scf.in` or `pw.in` to be rerun there without path
  editing.
- Recover validated, manually completed initial SCF, displacement, and
  `dynmat.x` calculations when their completion state was not yet committed.

### Changed

- Flatten external execution directories under `calculations/`, using names
  such as `init_scf_001`, `disp_0001_x_m_001`, and `dynmat_001`.
- Advance the dataset-state format so calculations using the previous nested
  layout are rejected explicitly instead of being mixed with the new layout.
- Force `disk_io='nowf'` in displaced QE inputs to reduce scratch-space use,
  while leaving the initial SCF policy unchanged.
- Use the QE-style `-inp` input option for FDVIB calculation and thermo
  commands and for FDVIB-launched `pw.x` and `dynmat.x` processes.
- Stream charge-density digest calculation and reuse the verified reference
  digest, avoiding whole-file memory use and repeated source-file hashing.

## [0.3.1] - 2026-06-24

### Fixed

- Copy `paw.txt` alongside the reference charge density when seeding displaced
  SCF calculations, fixing `Error in routine read_scf (1): Reading PAW becsum`
  for PAW pseudopotential workflows.

## [0.3.0] - 2026-06-24

### Added

- Add a resumable single-command calculation workflow through
  `fdvib -in fdvib.in`, with operating-system locking, immutable dataset
  fingerprints, append-only attempts, completion markers, result digests,
  and recovery of validated uncommitted results.
- Run and validate an unperturbed reference SCF before displacement jobs, then
  seed every displaced SCF from an independent copy of its converged charge
  density using FDVIB-managed `startingpot='file'`.
- Save the last converged reference-SCF total energy in Hartree as portable
  `electronic_structure.dat` metadata.
- Add optional automatic `dynmat.x` execution controlled by `run_dynmat` and
  `dynmat_command`.
- Add `fdvib shm RESULTS_DIR` for validated Shermo-compatible `.shm`
  export. Gas exports remove 3, 5, or 6 rigid modes according to molecular
  degrees of freedom, preserve separately calculated components of split
  degenerate vibrations, and retain remaining internal imaginary modes.
- Add regression coverage for reference-density seeding, restart and recovery,
  delayed `dynmat.x` execution, single-atom SHM output, nonlinear molecules,
  and linear CO2 mode selection.

### Changed

- Replace the separate prepare, run, and analyze commands with the resumable
  `fdvib -in fdvib.in` calculation workflow.
- Replace separate `mpi_command` and `pw_command` settings with one complete
  launcher command in `pw_command`.
- Rename configuration fields `workdir` to `outdir` and `output_prefix` to
  `prefix`; require an explicit `run_dynmat` logical setting.
- Keep Molden mode and thermochemistry analysis as separate commands; the
  latter now accepts an explicit `-in thermo.in` argument.
- Make configuration parsing strict: reject unknown, duplicate, malformed,
  or invalid fields instead of silently accepting ambiguous input.
- Keep `asr='no'` for both gas and local calculations and cross-check this,
  along with the interaction-block model, before SHM export.
- Update calculation layouts, examples, command documentation, recovery
  semantics, and the Shermo compatibility reference for the new workflow.

## [0.2.2] - 2026-06-22

### Fixed

- Parse only the total-force block from verbose `pw.x` output instead of
  allowing later decomposed force contributions to overwrite it.
- Refresh derived `forces.dat` files when preserving existing QE output, so
  previously misparsed force tables are repaired without rerunning `pw.x`.

## [0.2.1] - 2026-06-22

### Changed

- Use `asr='no'` for both local and gas calculations, and omit the unused
  `filmol` setting from generated `dynmat.in` files.
- Express gas pressure with `pressure_atm` and convert atmospheres to pascals
  using `1 atm = 101325 Pa`.
- Remove the absolute rigid-body frequency cutoff from gas RRHO calculations;
  report the largest removed frequency as a diagnostic instead.

## [0.2.0] - 2026-06-22

### Changed

- Update physical constants to the CODATA 2022 recommended values.
- Use `pw.x` as the executable command in both configuration examples.
- Split the monolithic implementation into common, QE I/O, execution,
  analysis, thermochemistry, and CLI translation units.
- Select gas-phase RRHO vibrations by molecular degrees of freedom: `3N-5`
  for linear molecules and `3N-6` for nonlinear molecules.
- Reject gas RRHO results with unresolved rigid-body modes, true imaginary
  vibrations, or a frequency-floor model.

### Added

- Report the rotor type and rigid-body/vibrational mode counts in gas
  thermochemistry output.
- Add end-to-end regression tests for gas RRHO mode selection and the
  prepare/run/restart workflow.

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

[Unreleased]: https://github.com/wxia529/fdvib/compare/v0.4.3...HEAD
[0.4.3]: https://github.com/wxia529/fdvib/compare/v0.4.2...v0.4.3
[0.4.2]: https://github.com/wxia529/fdvib/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/wxia529/fdvib/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/wxia529/fdvib/compare/v0.3.4...v0.4.0
[0.3.4]: https://github.com/wxia529/fdvib/compare/v0.3.3...v0.3.4
[0.3.3]: https://github.com/wxia529/fdvib/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/wxia529/fdvib/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/wxia529/fdvib/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/wxia529/fdvib/compare/v0.2.2...v0.3.0
[0.2.2]: https://github.com/wxia529/fdvib/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/wxia529/fdvib/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/wxia529/fdvib/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/wxia529/fdvib/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/wxia529/fdvib/releases/tag/v0.1.0
