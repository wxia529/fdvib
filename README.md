# FDVIB

FDVIB is a C++17 finite-difference vibration controller for
[Quantum ESPRESSO](https://www.quantum-espresso.org/). It generates Cartesian
displacements, runs `pw.x`, constructs a Gamma-point dynamical matrix from
central force differences, and can run QE's `dynmat.x`.

FDVIB supports:

- frozen-environment local harmonic calculations for selected atoms in a
  periodic system;
- isolated-molecule rigid-rotor harmonic-oscillator (RRHO) thermochemistry;
- compact Molden output for normal-mode visualization;
- Gaussian-like fake output for loading FDVIB modes in GaussView;
- Shermo 2.6.2-compatible `.shm` export;
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
thermo.in    optional thermochemistry settings
```

The QE input must use `ibrav=0` and explicit supported units on
`ATOMIC_POSITIONS` and `CELL_PARAMETERS`. It must also define
`calculation='scf'`, `tprnfor=.true.`, and an `outdir`. Do not set
`startingpot='file'`: FDVIB
runs one unperturbed reference SCF and seeds every displaced calculation from
its converged charge density. Both non-HDF5 QE `charge-density.dat` and HDF5
`charge-density.hdf5` outputs are supported.

Configuration templates are available in [`examples/local`](examples/local)
and [`examples/gas`](examples/gas). Gas calculations require
`selected_atoms = all`, an explicit spin multiplicity, and consistent QE spin
settings. Gas RRHO thermochemistry removes exactly five rigid-body modes for a
linear molecule or six for a nonlinear molecule and rejects any imaginary
frequency remaining in the true vibrational modes. Gas pressure is specified
with `pressure_atm`, where `1 atm = 101325 Pa`; gas calculations do not use
`zero_tolerance_cm1`.

## Workflow

From the calculation directory:

```sh
fdvib -inp fdvib.in
fdvib modes fdvib/results
fdvib fakeg fdvib/results
fdvib thermo fdvib/results -inp thermo.in
fdvib shm fdvib/results
```

The calculation command runs or resumes the reference SCF, all positive and
negative displacement SCFs, Hessian assembly, and optional `dynmat.x`
execution. Set `run_dynmat = true` or `false` in `fdvib.in`. Each displaced
SCF uses its own QE `outdir` and an FDVIB-injected `startingpot='file'`.
Displaced inputs also use `disk_io='nowf'` to avoid retaining unnecessary
wavefunction files.
Post-processing commands (`modes`, `fakeg`, `thermo`, and `shm`) may be rerun;
they overwrite only the files they generate and still reject damaged or
incomplete inputs.
Generated inputs use a calculation-local `outdir='./out'`. Runs are stored in
one flat directory level, for example `calculations/init_scf_001` and
`calculations/disp_0001_x_m_001`. An incomplete directory can be rerun there with
`pw.x -inp pw.in` (or `pw.x -inp scf.in` for the reference).
Completed directories are immutable snapshots and should be copied before
manual experiments.
Commands normally use `pw.x` from `PATH`; any file paths embedded in
`pw_command` should be absolute because the command runs inside that directory.
Generated local and gas inputs both use `asr='no'`; rigid-body modes are
excluded later according to molecular degrees of freedom.

Repeated calculation commands validate immutable dataset state and skip
completed stages. Failed attempts are retained instead of overwritten.
The stored electronic energy is the last converged `! total energy` from the
unperturbed reference SCF, converted from Ry to Hartree without adding ZPE.

Depending on the requested calculation and analysis commands, the result
directory contains:

```text
fdvib/results/
  system.dynG
  system.freq.out
  system.mold
  system_fake.out
  system.shm
  dynmat.in
  metadata.dat
  thermo.dat
```

The output prefix may differ from `system` according to `fdvib.in`.

## Documentation

See the [FDVIB reference](docs/index.md) for configuration fields, physical
models, output definitions, restart behavior, and diagnostics.

**Thermochemistry post-processing**: The mature and comprehensive
[Shermo](http://sobereva.com/soft/shermo) program is the recommended tool
for molecular thermochemistry analysis. FDVIB's `shm` command exports
results in Shermo's native format.
> *Quoted from http://sobereva.com/soft/shermo:* Shermo is a free, general,
> very easy-to-use and flexible code for calculating molecular thermochemistry
> data based on ideal gas assumption.

> If Shermo is utilized in your work, the following paper must be cited:
> Tian Lu, Qinxue Chen, *Shermo: A general code for calculating molecular
> thermodynamic properties*, Comput. Theor. Chem., 1200, 113249 (2021)
> DOI: [10.1016/j.comptc.2021.113249](https://doi.org/10.1016/j.comptc.2021.113249)

FDVIB's built-in `thermo` command is an alternative for quick checks.

Release history is recorded in the [changelog](CHANGELOG.md).

## License

FDVIB is distributed under the [BSD 3-Clause License](LICENSE).
