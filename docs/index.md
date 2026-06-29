# FDVIB reference

FDVIB is a tool for finite-difference vibrational analysis with Quantum
ESPRESSO. Start or resume a calculation with:

```bash
fdvib -inp fdvib.in
```

Normal-mode visualization and thermochemistry remain separate analyses:

```bash
fdvib modes fdvib/results
fdvib thermo fdvib/results -inp thermo.in
fdvib shm fdvib/results
```

For method details, see [Theory and method](theory.md).

## Calculation model

The calculation command performs these stages in order:

```text
validate inputs and dataset state
        |
unperturbed reference SCF
        |
copy the converged density into each displacement attempt
        |
positive and negative displacement SCFs with startingpot='file'
        |
central finite differences and Hessian symmetrization
        |
Gamma-point dynG and dynmat.in
        |
optional dynmat.x
```

There is no restart option to set. Run the same command again: FDVIB checks
the existing work and resumes from the first incomplete stage.

## Project files

A calculation starts with these user files:

```text
case_directory/
|-- scf.in          QE input for the original, unperturbed structure
|-- fdvib.in        FDVIB calculation settings
`-- thermo.in       optional, used only by the separate thermo command
```

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

Do not set this in `scf.in`:

```fortran
startingpot = 'file'
```

FDVIB rejects that value before starting a calculation. The unperturbed
reference SCF uses the starting-potential policy from `scf.in` (normally QE's
default atomic superposition). FDVIB injects `startingpot='file'` only into
the generated displaced inputs.

Every displaced calculation has an independent QE `outdir`. Before execution,
FDVIB copies the converged reference charge density into that calculation's
`<qe-prefix>.save` directory and verifies the copied data. Tasks never share a
writable QE scratch directory. FDVIB also inserts or replaces `disk_io='nowf'`
in displaced inputs to avoid storing wavefunction files that are not needed for
finite-difference forces. The initial SCF retains the user's `disk_io` policy.

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
| `selected_atoms` | One-based active atoms, or `all` for gas |
| `displacement_angstrom` | Positive displacement magnitude |
| `multiplicity` | Required positive spin multiplicity for gas |
| `pw_command` | Complete launcher and `pw.x` command |
| `prefix` | FDVIB result filename prefix |
| `run_dynmat` | Required logical controlling `dynmat.x` execution |
| `dynmat_command` | Complete `dynmat.x` command |

For gas calculations, the QE spin settings must agree with multiplicity. A
non-singlet requires `nspin=2` and
`tot_magnetization=multiplicity-1`.

## Initial SCF and displacement calculations

All external calculation runs use one flat directory level:

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
writes a nonempty charge-density file. For PAW calculations, the accompanying
`paw.txt` is also checked on restart.
FDVIB reads the last `! total energy = ... Ry` record from the unperturbed
reference SCF, divides it by two to convert Rydberg to Hartree, and stores the
pure reference value in result metadata. It does not use a displaced energy or
add ZPE. For calculations using electronic smearing, the QE total energy can
contain a smearing contribution; isolated-molecule inputs should avoid
artificial electronic-temperature contributions.

A failed or interrupted calculation remains available for diagnosis. The next
run creates a new numbered directory and starts again from the reference
density.

Generated `scf.in` and `pw.in` files are self-contained for diagnosis: FDVIB
starts QE inside each attempt directory, uses `outdir='./out'`, rewrites
relative `pseudo_dir` paths for that location, and redirects `wfcdir` to the
same isolated `./out`. Users can enter an attempt directory and run
`pw.x -inp pw.in` directly. Copy a completed directory before experimenting
with it. If an interrupted FDVIB process left behind a complete QE result,
FDVIB can check and reuse it on the next run.

Because `pw_command` is executed from the calculation directory, executables and
auxiliary files named in that command must either be discoverable through the
environment (for example `pw.x` on `PATH`) or use absolute paths.

In each `pw.out`, FDVIB reads the total-force table below the final
`Forces acting on atoms` heading. It ignores later tables that list individual
force contributions. A usable output must contain one total force for every
atom and `JOB DONE`, with no SCF non-convergence or QE error message.

## Hessian and QE dynamical matrix

For the active atom set, FDVIB evaluates the central finite difference:

```text
H[jβ,iα] = - (F[jβ](+δ[iα]) - F[jβ](-δ[iα])) / (2δ)
```

It reports the maximum antisymmetric component before symmetrization, then
writes:

```text
fdvib/results/<prefix>.dynG
fdvib/results/dynmat.in
fdvib/results/metadata.dat
```

Local calculations set `remove_interaction_blocks=.true.`. Gas calculations
retain the complete molecular Hessian. Both use `asr='no'`; rigid translations
and rotations are handled later during analysis.

With `run_dynmat=true`, FDVIB runs `dynmat_command` in an isolated attempt,
parses all `3N` frequencies and eigenvectors, and writes:

```text
fdvib/results/dynmat.out
fdvib/results/<prefix>.freq.out
```

With `run_dynmat=false`, calculation stops after `dynG` and `dynmat.in`.
Changing only this setting to `true` and repeating `fdvib -inp fdvib.in`
skips reference, displacement, and Hessian stages and runs only `dynmat.x`.

## Result metadata

FDVIB writes result-level metadata to:

```text
fdvib/results/metadata.dat
```

This file records post-processing information that is not contained in
`<prefix>.dynG` or `<prefix>.freq.out`:

```text
program = qe
mode_selection = gas
electronic_energy_hartree = -76.0
multiplicity = 1
selected_atoms = all
```

`program` is reserved for source-program dispatch; the current reader accepts
only `qe`. `mode_selection` controls how exporters select frequencies:
`all` writes all nonzero frequencies, `gas` removes rigid translations and
rotations by molecular geometry, and `local` keeps `3 * selected_atoms` modes
for SHM export. `electronic_energy_hartree` is the unperturbed SCF energy used
by `.shm`; if the field is absent, SHM writes `0.0`. `multiplicity` is written
as the ground-state electronic degeneracy in `.shm`. `selected_atoms` records
the one-based active atom list or `all`.

The file is generated by the calculation workflow. Users may edit it before
running post-processing commands when they intentionally want to change export
metadata, for example `mode_selection`, `multiplicity`, or the electronic
energy used in `.shm`. Editing `metadata.dat` does not change the completed
force calculation or the Hessian; it only affects later exports. If
`metadata.dat` is missing, exporters use these defaults:
`program=qe`, `mode_selection=all`, `multiplicity=1`, and
`electronic_energy_hartree=0.0`.

## State and recovery

To continue an interrupted calculation, run the same `fdvib -inp fdvib.in`
command again with the same `outdir`. FDVIB checks that existing results belong
to the same calculation before reusing them.

You may change execution-only settings such as `pw_command`, `dynmat_command`,
and `run_dynmat`. For example, a calculation first run with
`run_dynmat=false` can later be rerun with `run_dynmat=true`; completed SCF and
Hessian stages are reused and only `dynmat.x` is added.

Do not reuse the same `outdir` after changing `scf.in`, `system_type`,
`selected_atoms`, `displacement_angstrom`, `multiplicity`, or the FDVIB result
`prefix`. Choose a new `outdir` for that new calculation. This prevents results
from two different physical setups from being mixed. The internal
`fdvib/state/dataset.state` file records the original calculation settings for
this check and should not be edited.

FDVIB prevents two processes from writing to the same `outdir` at the same
time. Before reusing a completed stage, it checks that the recorded files are
still present and unchanged. Restart behavior is summarized below:

| State | Behavior |
|---|---|
| No attempt | Run it |
| Failed or interrupted attempt | Preserve it and create a new attempt |
| Valid completed stage | Skip it |
| Complete marker with missing or changed result | Stop and report corruption |
| Fully written result without final marker | Validate and recover the marker |
| Partial uncommitted publication | Preserve it under `fdvib/failed` and retry |

## Normal-mode analysis

Run after a calculation with valid `dynmat.x` results:

```bash
fdvib modes fdvib/results
```

The command reads `<prefix>.dynG`, `<prefix>.freq.out`, and the optional
`metadata.dat`, then writes a CP2K-style `<prefix>.mol`.
It preserves imaginary-frequency signs, includes the complete geometry, and
writes zero displacement for frozen atoms. Re-running the command overwrites
the generated `.mol` file.

For `mode_selection=gas`, FDVIB writes only the molecular internal modes. For
`local` and `all`, only frequencies stored as exactly `0.0` are omitted; small
nonzero frequencies are retained. Every export includes the three QE cell
vectors in the Molden `[Cell]` section;
`gas` controls mode selection, not whether the underlying periodic supercell
exists. Atom and `FR-COORD` coordinates are written in Bohr, while `[Cell]`
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
It also reads `metadata.dat` when present; see the result metadata section for
the recorded fields and defaults. The default `all` mode writes every nonzero
frequency and omits only frequencies exactly equal to `0.0`.

The command writes `<prefix>.shm` with electronic energy, wavenumbers, atoms,
and electronic levels. Coordinates are written in Angstrom, masses in amu, and
the ground-state electronic degeneracy is the configured multiplicity.

QE species labels such as `C1` or `O_ads` are reduced to their leading
standard element symbol; labels that cannot be mapped to H--Og are rejected.

For gas systems, FDVIB uses the same inertia threshold as Shermo to classify
an atom, linear molecule, or nonlinear molecule. It removes the 3, 5, or 6
frequencies closest to zero and writes the remaining internal vibrations:

```text
atom       0 modes
linear     3N-5 modes
nonlinear  3N-6 modes
```

Numerically split degenerate modes are preserved as separate original
frequencies. For example, the two CO2 bending components are both written;
FDVIB does not average them. A true transition-state imaginary frequency is
also retained when it is not among the rigid modes closest to zero.

Mode identification is frequency-based. The export summary reports both the
largest removed absolute frequency and the smallest retained absolute
frequency. If these groups are not clearly separated, inspect the modes and
improve geometry/force convergence before using the thermochemistry; a very
soft true vibration can otherwise be confused with a residual rigid mode.

For local systems, FDVIB discards the `3(N-N_active)` modes closest to zero,
writes the remaining `3N_active` modes and all atoms. Run Shermo as:

```bash
Shermo fdvib/results/system.shm -imode 1 -PGlabel C1
```

For gas `.shm` files, run Shermo directly according to your target
analysis and the Shermo documentation, for example:

```bash
Shermo fdvib/results/molecule.shm
```

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

Local analysis reports vibrational ZPE, internal energy, entropy, and free
energy. Gas RRHO additionally includes translation, rotation, and electronic
degeneracy. It removes three rigid modes for an atom, five for a linear
molecule, or six for a nonlinear molecule, then rejects any remaining
non-positive vibrational frequency. Gas RRHO requires `metadata.dat` with
`mode_selection = gas` so the spin multiplicity is unambiguous. Results are
written to `fdvib/results/thermo.dat`. The file reports each energy table in eV,
kcal/mol, and kJ/mol; the corresponding entropy units are eV/K,
kcal/(mol K), and kJ/(mol K).

Changing temperature or other thermochemistry settings requires only another
`fdvib thermo` command; it never reruns the electronic-structure calculation.

## Acknowledgements

FDVIB uses [Quantum ESPRESSO](https://www.quantum-espresso.org/) as its
electronic-structure backend, and its `.shm` exporter follows the documented
input format of [Shermo](http://sobereva.com/soft/shermo/). We acknowledge the
developers and contributors of both projects for their work.
