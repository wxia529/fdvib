# FDVIB

FDVIB is a command-line tool written in Go, designed for local
finite-difference vibrational analysis in periodic systems with
[Quantum ESPRESSO](https://www.quantum-espresso.org/). It can also treat
isolated molecules. FDVIB generates Cartesian displacements, runs `pw.x`,
constructs a Gamma-point dynamical matrix from central force differences, and
can run QE's `dynmat.x`.

FDVIB supports:

- local harmonic calculations for selected atoms in a periodic system;
- isolated-molecule rigid-rotor harmonic-oscillator (RRHO) thermochemistry;
- Molden-style output for normal-mode visualization;
- Shermo-compatible `.shm` export;
- harmonic treatment of positive frequencies, with an option to raise low
  frequencies to a configurable minimum.

It is not a replacement for `ph.x`, `q2r.x`, Phonopy, or full periodic phonon
calculations.

## Requirements

- Go 1.22 or later
- Quantum ESPRESSO `pw.x` and `dynmat.x`

## Diagnostics

To diagnose slow task transitions or failures, add `--debug` to a command:

```sh
fdvib --debug -inp fdvib.in
```

Debug progress goes to stderr; a unique `fdvib-debug-*.jsonl` file is created
in the current directory. Use `--debug-log FILE` with `--debug` to choose a
new log file. See [diagnostics](docs/index.md#diagnostics) for timing details.

## Build

```sh
go build -ldflags "-X main.version=1.1.2" -o fdvib ./cmd/fdvib
```

The executable is written to `./fdvib` (or `bin/fdvib` with
`go build ./cmd/fdvib`). To install it:

```sh
go install -ldflags "-X main.version=1.1.2" github.com/wxia529/fdvib/cmd/fdvib
```

## Linux releases

Version tags produce a prebuilt Linux x86_64 archive on the GitHub Releases
page. The binary is a statically linked Go executable with no glibc
dependency; it runs on any modern Linux distribution. The archive also
contains the reference documentation, license, and configuration examples.
Quantum ESPRESSO is not bundled.

```sh
tar -xzf fdvib-X.Y.Z-linux-x86_64.tar.gz
install fdvib-X.Y.Z-linux-x86_64/bin/fdvib "$HOME/.local/bin/fdvib"
```

## Input

A typical calculation directory contains:

```text
scf.in       Quantum ESPRESSO SCF input
fdvib.in     displacement and execution settings
thermo.in    optional thermochemistry settings
```

The QE input does not have to be named `scf.in`; set `scf_input` in
`fdvib.in` to the actual filename.

Create a starter `fdvib.in` in the current directory with:

```sh
fdvib init local
# or
fdvib init gas
```

The command creates only `fdvib.in` and refuses to overwrite an existing
file. Review the generated settings before starting a calculation.

The QE input must use `ibrav=0` and explicit supported units on
`ATOMIC_POSITIONS` and `CELL_PARAMETERS`. It must also define
`calculation='scf'`, `tprnfor=.true.`, and an `outdir`. Do not set
`startingpot='file'`: FDVIB runs one unperturbed reference SCF and seeds every
displaced calculation from its converged charge density. FDVIB detects the
charge-density representation written by the installed QE build.

Configuration templates are available in [`examples/local`](examples/local)
and [`examples/gas`](examples/gas). Gas calculations require
`selected_atoms = all`, an explicit spin multiplicity, and consistent QE spin
settings. The [reference documentation](docs/index.md) describes the supported
input fields and thermochemistry models.

## Workflow

From the calculation directory:

```sh
fdvib -inp fdvib.in
fdvib modes fdvib/results
fdvib thermo fdvib/results -inp thermo.in
fdvib shm fdvib/results
```

The first command runs the reference SCF, the positive and negative
displacements, force-constant matrix assembly, and optionally `dynmat.x`.
Running it again with the same input resumes from the first incomplete stage.
Failed attempts are kept for diagnosis.

The other three commands are independent analyses. They may be run repeatedly
and replace only their own output files. See the reference documentation for
restart rules, generated QE directories, and manual diagnosis.

The commands above produce some or all of these files:

```text
fdvib/results/
  system.dynG
  system.freq.out
  system.mol
  system.shm
  dynmat.in
  metadata.dat
  thermo.dat
```

The output prefix may differ from `system` according to `fdvib.in`.

## Documentation

See the [FDVIB reference](docs/index.md) for configuration fields, output
definitions, restart behavior, and diagnostics. The physical method is
described in [Theory and method](docs/theory.md).

The `shm` command writes input for Shermo.

**If Shermo is utilized in your work, the following paper must be cited:**

Tian Lu, Qinxue Chen, *Shermo: A general code for calculating molecular
thermodynamic properties*, Comput. Theor. Chem., 1200, 113249 (2021). DOI:
[10.1016/j.comptc.2021.113249](https://doi.org/10.1016/j.comptc.2021.113249)

For usage instructions and further information, see the
[Shermo official website](http://sobereva.com/soft/shermo/).

Release history is recorded in the [changelog](CHANGELOG.md).

## Acknowledgements

FDVIB uses [Quantum ESPRESSO](https://www.quantum-espresso.org/) as its
electronic-structure backend, and its `.shm` exporter follows the documented
input format of [Shermo](http://sobereva.com/soft/shermo/). We acknowledge the
developers and contributors of both projects for their work.

## License

Copyright (c) 2026 Wanting Xia. FDVIB is distributed under the
[BSD 3-Clause License](LICENSE).
