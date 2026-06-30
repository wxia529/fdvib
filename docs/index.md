# FDVIB reference

FDVIB is designed for local vibrational analysis in periodic systems using
finite differences of Quantum ESPRESSO forces. It displaces selected atoms
while keeping the surrounding atoms as an unselected environment. It can
also calculate isolated-molecule vibrations by displacing every atom.

FDVIB is not a replacement for a full periodic phonon calculation with
`ph.x`. For the equations and mode-selection rules used here, see
[Theory and method](theory.md).

## Calculation types

`local` is for periodic systems or embedded models. FDVIB displaces only the
selected atoms, which form the active region specified by `selected_atoms`.
It constructs a local harmonic approximation for that region. The full
structure is retained for visualization and export, but the result is not the
complete vibrational spectrum of the periodic system.

`gas` is for isolated molecules. FDVIB displaces every atom and constructs the
full molecular force-constant matrix. Translation and rotation are removed
during gas-mode post-processing.

## Quick start

Start or resume a calculation with:

```bash
fdvib -inp fdvib.in
```

Run analyses separately after `dynmat.x` results are available:

| Purpose | Command | Result |
|---|---|---|
| Visualize modes | `fdvib modes fdvib/results` | Molden-style `.mol` file |
| Use Shermo | `fdvib shm fdvib/results` | `.shm` file for Shermo |
| Compute thermochemistry | `fdvib thermo fdvib/results -inp thermo.in` | `thermo.dat` |

A typical workflow is:

1. prepare and test the QE SCF input with `pw.x`;
2. run `fdvib init local` or `fdvib init gas`, then edit the generated
   `fdvib.in`;
3. run `fdvib -inp fdvib.in`;
4. run `fdvib modes`, `fdvib shm`, or `fdvib thermo` as needed.

The calculation command checks the input, runs the reference SCF and all
positive and negative displacement calculations, builds the force-constant
matrix, and optionally runs `dynmat.x`. There is no restart option to set. Run
the same command again to continue from the first incomplete stage.

## Required input files

A calculation starts with these user files:

```text
case_directory/
|-- scf.in          QE input for the original, unperturbed structure
|-- fdvib.in        FDVIB calculation settings
`-- thermo.in       optional, used only by the separate thermo command
```

The tree uses `scf.in` as an example. The QE input may have any filename as
long as `scf_input` in `fdvib.in` points to it.

To create a starter `fdvib.in` in the current directory, run one of:

```bash
fdvib init local
fdvib init gas
```

The command creates only `fdvib.in`; it does not create `scf.in` or
`thermo.in`. It refuses to overwrite an existing `fdvib.in`. The generated
`scf_input = scf.in` is only a default; change it if the QE input has another
name.

In `fdvib.in`, `scf_input = scf.in` means the `scf.in` beside `fdvib.in`:

```text
scf_input = scf.in
outdir = fdvib
```

Run FDVIB from `case_directory`:

```bash
fdvib -inp fdvib.in
```

FDVIB leaves the top-level `case_directory/scf.in` unchanged. It creates new QE
inputs for execution under `fdvib/calculations/`:

```text
fdvib/calculations/
|-- init_scf_001/scf.in       generated input, no displacement
|-- disp_0001_x_p_001/pw.in  generated input, atom 1 displaced along +x
`-- disp_0001_x_m_001/pw.in  generated input, atom 1 displaced along -x
```

The generated `pw.in` files contain the displaced coordinates; the user's
top-level `scf.in` does not.

## QE input requirements

`scf.in` must define a fixed-ion SCF calculation with:

```fortran
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  prefix = 'pwscf'
  outdir = './out'
/
&SYSTEM
  ibrav = 0
  nat = ...
  ntyp = ...
/

ATOMIC_POSITIONS angstrom
...

CELL_PARAMETERS angstrom
...
```

FDVIB requires `ibrav=0` and explicit units on the position and cell cards.
Accepted position units are `angstrom`, `bohr`, `alat`, and `crystal`.
Accepted cell units are `angstrom`, `bohr`, and `alat`. QE spellings with
bare units, parentheses, or braces are accepted, for example
`ATOMIC_POSITIONS angstrom`, `ATOMIC_POSITIONS (bohr)`, and
`ATOMIC_POSITIONS {crystal}`. `crystal_sg` and omitted/deprecated default units
are not supported. `tprnfor=.true.` is required for force extraction.

FDVIB first creates a reference input from the user's original `scf.in`,
without displacing any atoms. For each displaced calculation, it then copies
the converged reference density and inserts `startingpot='file'`. This gives
the positive and negative displacements a consistent starting point and
reduces unnecessary SCF work.

Therefore, do not set the following value in the original `scf.in`:

```fortran
startingpot = 'file'
```

FDVIB rejects it before starting. Each displaced calculation has an
independent QE `outdir`, so calculations never share writable QE scratch data.
FDVIB also inserts or replaces `disk_io='nowf'` because wavefunction files are
not required for the finite-difference forces. The reference SCF retains the
user's `disk_io` policy.

The original `scf.in` should be tested independently. Tight electronic
convergence is important because force noise directly affects frequencies.

## fdvib.in

Local periodic example:

```text
scf_input = scf.in
outdir = fdvib
system_type = local
selected_atoms = 65,66,67
displacement_angstrom = 0.01
pw_command = mpirun -np 8 pw.x
prefix = system
run_dynmat = true
dynmat_command = dynmat.x
```

Gas molecule example:

```text
scf_input = scf.in
outdir = fdvib
system_type = gas
selected_atoms = all
displacement_angstrom = 0.01
multiplicity = 1
pw_command = mpirun -np 8 pw.x
prefix = molecule
run_dynmat = true
dynmat_command = dynmat.x
```

`fdvib.in` is a plain key-value file, not a QE-style namelist. Write one
`key = value` assignment per line. The parser ignores `#` and `!` comments and
rejects unknown or duplicate parameters. Values may be quoted, but quotes are
not required.

| Parameter | Meaning |
|---|---|
| `scf_input` | External QE SCF input, relative to `fdvib.in` |
| `outdir` | FDVIB calculation directory |
| `system_type` | `local` or `gas` |
| `selected_atoms` | One-based selected atoms, or `all` for gas |
| `displacement_angstrom` | Positive displacement magnitude |
| `multiplicity` | Required positive spin multiplicity for gas |
| `pw_command` | Complete launcher and `pw.x` command |
| `prefix` | FDVIB result filename prefix |
| `run_dynmat` | Whether to run `dynmat.x` after building the force-constant matrix |
| `dynmat_command` | Complete `dynmat.x` command |

The FDVIB `prefix` is not the QE `prefix`. QE's `prefix` names the
`<qe-prefix>.save` directory used by `pw.x`. The FDVIB `prefix` names result
files such as `system.dynG`, `system.freq.out`, and `system.shm`. The two values
do not need to match.

For gas calculations, the QE spin settings must agree with multiplicity. A
non-singlet requires `nspin=2` and
`tot_magnetization=multiplicity-1`.

## What happens during a calculation

The calculation proceeds in five stages:

1. run one reference SCF at the input geometry;
2. run positive and negative displacements for every selected Cartesian
   coordinate;
3. obtain the Cartesian force-constant matrix from the force differences;
4. write a QE-compatible Gamma-point dynamical-matrix file and `dynmat.in`;
5. run `dynmat.x` when `run_dynmat=true`.

For each selected displacement, the positive and negative force calculations
provide one column of the force-constant matrix. The exact matrix definition
and sign convention are given in [Theory and method](theory.md).

Each QE attempt has its own directory directly under `fdvib/calculations/`:

```text
fdvib/calculations/
|-- init_scf_001/
|   |-- scf.in
|   |-- scf.out
|   `-- out/<qe-prefix>.save/
|       |-- charge-density.*
|       `-- paw.txt                 # present for PAW calculations
|-- disp_0065_x_p_001/
|   |-- pw.in
|   |-- pw.out
|   |-- forces.dat
|   `-- out/<qe-prefix>.save/
|       |-- charge-density.*
|       `-- paw.txt                 # copied for PAW calculations
`-- disp_0065_x_m_001/
```

`init_scf` is the SCF at the input geometry. A displacement name records the
one-based QE atom index, Cartesian direction, and positive (`p`) or negative
(`m`) displacement. The final three digits are the execution number; a failed
run is retained and the retry receives `_002`. FDVIB records which numbered
directory completed successfully and checks its files before reusing it.

The initial SCF is accepted only after QE exits successfully, reports
`JOB DONE`, reaches SCF convergence, prints one total force for every atom, and
writes a nonempty charge-density file. For PAW calculations, FDVIB also checks
`paw.txt` on restart.

FDVIB stores the last `! total energy = ... Ry` value from the unperturbed SCF
and converts it from Rydberg to Hartree. Displaced energies and ZPE are not
included. This is expected: FDVIB obtains force constants from force
differences, not from displaced total-energy differences. ZPE and thermal
corrections are computed later by `fdvib thermo` or by Shermo. If the QE
calculation uses electronic smearing, the reported total energy may include a
smearing contribution. This is usually undesirable for an isolated molecule.

A failed or interrupted calculation remains available for diagnosis. The next
run creates a new numbered directory and starts again from the reference
density.

Generated `scf.in` and `pw.in` files can be rerun from their attempt
directories. FDVIB makes this possible by:

- starting QE inside the attempt directory;
- setting `outdir='./out'`;
- rewriting a relative `pseudo_dir` for the new location;
- redirecting `wfcdir` to the same `./out` directory.

Enter an attempt directory and run `pw.x -inp pw.in` to diagnose a displaced
calculation. Use `pw.x -inp scf.in` for the initial SCF. Copy a completed
directory before experimenting with it. If FDVIB was interrupted after QE
finished, the next run can check and reuse that result.

`pw_command` also runs inside the attempt directory. Programs such as `pw.x`
must therefore be available on `PATH`, and other file paths in the command
should be absolute.

In each `pw.out`, FDVIB reads the total-force table below the final
`Forces acting on atoms` heading. It ignores later tables that list individual
force contributions. A usable output must contain one total force for every
atom and `JOB DONE`, with no SCF non-convergence or QE error message.

## Main calculation outputs

FDVIB builds the force-constant matrix from the positive and negative force
calculations. It reports the maximum antisymmetric component before
symmetrizing the matrix, which helps identify noisy or poorly converged
forces.

The calculation then writes:

```text
fdvib/results/<prefix>.dynG
fdvib/results/dynmat.in
fdvib/results/metadata.dat
```

`<prefix>.dynG` is generated by FDVIB from the finite-difference force
constants. It uses a QE-compatible Gamma-point dynamical-matrix format, but it
is not output from `ph.x`.

Local calculations set `remove_interaction_blocks=.true.`; see
[Theory and method](theory.md#qe-dynamical-matrix-file) for how the zero-padded
local matrix is handled during QE post-processing. Gas calculations retain the
complete molecular force-constant matrix. Both use `asr='no'`; rigid
translations and rotations are handled later during analysis.

With `run_dynmat=true`, FDVIB runs `dynmat_command`, parses all `3N` frequencies
and eigenvectors, and writes:

```text
fdvib/results/dynmat.out
fdvib/results/<prefix>.freq.out
```

With `run_dynmat=false`, calculation stops after `dynG` and `dynmat.in`.
Changing only this setting to `true` and repeating `fdvib -inp fdvib.in`
skips reference, displacement, and force-constant stages and runs only
`dynmat.x`.

## Normal-mode analysis

Run after a calculation with valid `dynmat.x` results:

```bash
fdvib modes fdvib/results
```

The command reads `<prefix>.dynG`, `<prefix>.freq.out`, and the optional
`metadata.dat`, then writes a Molden-style `<prefix>.mol` using the same
coordinate and cell layout convention as CP2K.
It preserves imaginary-frequency signs, includes the complete geometry, and
writes zero displacement for unselected environment atoms. Re-running the
command overwrites the generated `.mol` file.

For `mode_selection=gas`, FDVIB writes only the molecular internal modes. For
`local` and `all`, only frequencies stored as exactly `0.0` are omitted; small
nonzero frequencies are retained. Every export includes the three QE cell
vectors in the Molden `[Cell]` section. The `gas` setting changes only the
selected modes; it does not remove the QE supercell. Atom and `FR-COORD`
coordinates are written in Bohr, while `[Cell]`
vectors are written in Angstrom, matching the CP2K layout. The `[INT]` section
is omitted because FDVIB does not calculate IR intensities.

QE's `dynmat.x` can also write `dynmat.mold` through its `filmol` option.
FDVIB writes its own `.mol` from the same normalized modes so that signed
imaginary frequencies, cell and atom information, and the selected set of
modes are preserved. It does not diagonalize the dynamical matrix again.

## Shermo `.shm` export

After `dynmat.x` has produced a frequency file, run:

```bash
fdvib shm fdvib/results
```

The command produces an input file for Shermo.

**If Shermo is utilized in your work, the following paper must be cited:**

Tian Lu, Qinxue Chen, *Shermo: A general code for calculating molecular
thermodynamic properties*, Comput. Theor. Chem., 1200, 113249 (2021). DOI:
[10.1016/j.comptc.2021.113249](https://doi.org/10.1016/j.comptc.2021.113249)

For usage instructions and further information, see the
[Shermo official website](http://sobereva.com/soft/shermo/).

The command requires exactly one `<prefix>.dynG` and one `<prefix>.freq.out`.
It also reads `metadata.dat` when present; see
[Result metadata](#result-metadata) for the recorded fields and defaults. The
default `all` mode writes every nonzero frequency and omits only frequencies
exactly equal to `0.0`.

The command writes `<prefix>.shm` with electronic energy, wavenumbers, atoms,
and electronic levels. Coordinates are written in Angstrom, masses in amu, and
the ground-state electronic degeneracy is the configured multiplicity.

QE species labels such as `C1` or `O_ads` are reduced to their leading
standard element symbol. Labels that cannot be mapped to an element from H
through Og are rejected.

For gas systems, FDVIB uses the principal moments of inertia to classify an
atom, linear molecule, or nonlinear molecule. It removes the 3, 5, or 6
frequencies closest to zero, respectively, and writes the remaining internal
vibrations. The detailed selection rule is given in
[Theory and method](theory.md#gas-calculations).

Numerically split degenerate modes are preserved as separate original
frequencies. For example, the two CO2 bending components are both written;
FDVIB does not average them. A true transition-state imaginary frequency is
also retained when it is not among the rigid modes closest to zero.

Mode identification is frequency-based. The export summary reports both the
largest removed absolute frequency and the smallest retained absolute
frequency. If these groups are not clearly separated, inspect the modes and
improve geometry/force convergence before using the thermochemistry; a very
soft true vibration can otherwise be confused with a residual rigid mode.

For a local system with $N_{\mathrm{active}}$ selected atoms, FDVIB removes
the modes closest to zero until $3N_{\mathrm{active}}$ remain, then writes
those modes and all atoms to the SHM file. Run Shermo as:

```bash
Shermo fdvib/results/system.shm -imode 1 -PGlabel C1
```

Here `-imode 1` is used for the local vibrational treatment, while
`-PGlabel C1` avoids assigning molecular symmetry to the periodic structure.

The resulting Shermo correction is a local vibrational correction for the
selected active region, not the thermochemistry of the full periodic
structure.

For gas `.shm` files, Shermo can usually be run directly:

```bash
Shermo fdvib/results/molecule.shm
```

Additional Shermo options are described on its official website.

Re-running the command safely replaces the previously generated `.shm` file.

## Built-in thermochemistry

Run the built-in thermo analysis with:

```bash
fdvib thermo fdvib/results -inp thermo.in
```

Local harmonic example:

```text
model = local_harmonic
temperature_k = 298.15
low_frequency_model = frequency_floor
frequency_floor_cm1 = 100.0
zero_tolerance_cm1 = 1.0
```

Gas RRHO example:

```text
model = gas_rrho
temperature_k = 298.15
pressure_atm = 1.0
symmetry_number = 1
electronic_degeneracy = auto
rotor_type = auto
low_frequency_model = harmonic
```

`local_harmonic` reports ZPE, internal energy, entropy, and free energy from
the selected local vibrational modes only. It does not include translation or
rotation of the full periodic system. `gas_rrho` uses the rigid-rotor
harmonic-oscillator model and includes translational, rotational, vibrational,
and electronic-degeneracy contributions. It removes three rigid modes for an
atom, five for a linear molecule, or six for a nonlinear molecule, then
rejects any remaining non-positive vibrational frequency. Gas RRHO requires
`metadata.dat` with
`mode_selection = gas`. The same file supplies the multiplicity used when
`electronic_degeneracy = auto`. Results are written to
`fdvib/results/thermo.dat`. The file reports each energy table in eV,
kcal/mol, and kJ/mol; the corresponding entropy units are eV/K,
kcal/(mol K), and kJ/(mol K).

Changing temperature or other thermochemistry settings requires only another
`fdvib thermo` command; it never reruns the electronic-structure calculation.

## State and recovery

To continue an interrupted calculation, run the same `fdvib -inp fdvib.in`
command again with the same `outdir`. FDVIB checks that existing results belong
to the same calculation before reusing them.

You may change execution-only settings such as `pw_command`, `dynmat_command`,
and `run_dynmat`. For example, a calculation first run with
`run_dynmat=false` can later be rerun with `run_dynmat=true`; completed SCF and
force-constant stages are reused and only `dynmat.x` is added.

Do not reuse the same `outdir` after changing `scf.in`, `system_type`,
`selected_atoms`, `displacement_angstrom`, `multiplicity`, or the FDVIB result
`prefix`. Choose a new `outdir` so results from different physical setups are
not mixed. The internal `fdvib/state/dataset.state` records the original
settings for this check and should not be edited.

FDVIB prevents two processes from writing to the same `outdir` at the same
time. Before reusing a completed stage, it checks that the recorded files are
still present and unchanged.

| Situation | What FDVIB does |
|---|---|
| The stage has not started | Run it |
| A previous attempt failed or was interrupted | Keep that directory and start a numbered retry |
| The recorded result is complete and unchanged | Skip the stage |
| A previously completed file is missing or changed | Stop and report the affected file |
| QE finished before FDVIB was interrupted | Check the files and reuse the result if valid |
| Incomplete files were left in `results/` | Move them to `fdvib/failed` and rebuild them |

## Result metadata

FDVIB writes post-processing settings to:

```text
fdvib/results/metadata.dat
```

The file supplements information in `<prefix>.dynG` and
`<prefix>.freq.out`:

```text
program = qe
mode_selection = gas
electronic_energy_hartree = -76.0
multiplicity = 1
selected_atoms = all
```

| Field | Meaning |
|---|---|
| `program` | Electronic-structure program; currently only `qe` is supported |
| `mode_selection` | Selects modes for later analysis; the rules are described in [Theory and method](theory.md#export-mode-selection) |
| `electronic_energy_hartree` | Unperturbed SCF energy written to `.shm`; SHM uses `0.0` if it is absent |
| `multiplicity` | Ground-state electronic degeneracy written to `.shm` |
| `selected_atoms` | One-based selected atom list, or `all` |

The `program` field is reserved for future support of other
electronic-structure programs.

The calculation command generates this file. You may edit `mode_selection`,
`multiplicity`, or `electronic_energy_hartree` before running an analysis
command. These changes affect post-processing only; they do not alter the
completed force calculation or the force-constant matrix. If `metadata.dat`
is missing, the defaults are `program=qe`, `mode_selection=all`,
`multiplicity=1`, and `electronic_energy_hartree=0.0`.

## Acknowledgements

FDVIB uses [Quantum ESPRESSO](https://www.quantum-espresso.org/) as its
electronic-structure backend, and its `.shm` exporter follows the documented
input format of [Shermo](http://sobereva.com/soft/shermo/). We acknowledge the
developers and contributors of both projects for their work.
