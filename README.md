# FDVIB

FDVIB is a C++17 finite-difference vibration controller for
[Quantum ESPRESSO](https://www.quantum-espresso.org/). It generates Cartesian
displacements, runs `pw.x`, constructs a Gamma-point dynamical matrix from
central force differences, and prepares input for QE's `dynmat.x`.

FDVIB supports:

- frozen-environment local harmonic calculations for selected atoms in a
  periodic system;
- isolated-molecule rigid-rotor harmonic-oscillator (RRHO) thermochemistry;
- compact Molden output for normal-mode visualization;
- harmonic and frequency-floor treatments of low positive frequencies.

It is not a replacement for `ph.x`, `q2r.x`, Phonopy, or full periodic phonon
calculations.

## Requirements

- CMake 3.12 or later
- A C++17 compiler
- Quantum ESPRESSO `pw.x` and `dynmat.x`

## Build

```sh
cmake -S . -B build
cmake --build build -j
```

The executable is written to `build/fdvib`. To install it:

```sh
cmake --install build --prefix "$HOME/.local"
```

## Linux releases

Version tags produce a prebuilt Linux x86_64 archive on the GitHub Releases
page. The binary targets glibc 2.17 or newer and statically links `libstdc++`
and `libgcc`; it is therefore suitable for CentOS/RHEL 7 and most newer glibc
distributions. The archive also contains the reference documentation, license,
and configuration examples. Quantum ESPRESSO is not bundled.

```sh
tar -xzf fdvib-X.Y.Z-linux-x86_64-glibc-2.17.tar.gz
install fdvib-X.Y.Z-linux-x86_64-glibc-2.17/bin/fdvib "$HOME/.local/bin/fdvib"
```

## Input

A calculation directory contains:

```text
scf.in       Quantum ESPRESSO SCF input
fdvib.in     displacement and execution settings
thermo.in    thermochemistry settings
```

The QE input must use `ibrav=0`, `ATOMIC_POSITIONS angstrom`, and
`CELL_PARAMETERS angstrom`. It must also define `calculation='scf'`,
`tprnfor=.true.`, and an `outdir`.

Configuration templates are available in [`examples/local`](examples/local)
and [`examples/gas`](examples/gas). Gas calculations require
`selected_atoms='all'`, an explicit spin multiplicity, and consistent QE spin
settings. Gas RRHO thermochemistry removes exactly five rigid-body modes for a
linear molecule or six for a nonlinear molecule and rejects any imaginary
frequency remaining in the true vibrational modes.

## Workflow

From the calculation directory:

```sh
fdvib prepare fdvib.in
fdvib run fdvib.in
fdvib analyze fdvib.in

cd fdvib/results
dynmat.x -in dynmat.in > dynmat.out
cd ../..

fdvib modes fdvib/results
fdvib thermo fdvib/results
```

`prepare` creates positive and negative displacements for each selected
Cartesian coordinate. `run` executes the corresponding `pw.x` calculations.
`analyze` forms and symmetrizes the finite-difference Hessian and writes the QE
dynamical matrix and `dynmat.in`. Execution of `dynmat.x` remains an explicit
QE post-processing step.

Existing QE output and Molden files are not overwritten. Dataset settings used
by `analyze` are loaded from immutable snapshots created by `prepare`.

The result directory contains:

```text
fdvib/results/
  system.dynG
  system.freq.out
  system.mold
  dynmat.in
  fdvib.in.reference
  thermo.in
  thermo.dat
```

The output prefix may differ from `system` according to `fdvib.in`.

## Documentation

See the [FDVIB reference](docs/index.md) for configuration fields, physical
models, output definitions, restart behavior, and diagnostics.
Release history is recorded in the [changelog](CHANGELOG.md).

## License

FDVIB is distributed under the [BSD 3-Clause License](LICENSE).
