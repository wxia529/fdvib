# FDVIB reference

FDVIB is a finite-difference vibration controller for Quantum ESPRESSO. A
calculation runs from one command and can be resumed safely:

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

There is no user-selectable restart mode. Repeating the same command validates
completed stages and resumes the first incomplete stage. FDVIB never treats a
file as complete merely because it exists.

## Project files

A typical user working directory starts with:

```text
case_directory/
|-- scf.in
|-- fdvib.in
`-- thermo.in       optional, used only by the separate thermo command
```

`scf.in` stays external; it is not embedded in `fdvib.in`. During execution
FDVIB creates its own `fdvib/calculations/` subdirectory for generated QE
attempts; that internal name is separate from the user's working directory.

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
run is retained and the retry receives `_002`. Completion snapshots such as
`state/init_scf.complete` and `state/disp_0065_x_p.complete` identify the
accepted calculation directory and record its data digests.

The initial SCF is committed only after FDVIB verifies the process exit
status, `JOB DONE`, SCF convergence, the complete force block, and a nonempty
charge-density file. When QE produces `paw.txt`, FDVIB includes its presence
and digest in the completion snapshot and verifies both files on restart.
FDVIB reads the last `! total energy = ... Ry` record from the unperturbed
reference SCF, divides it by two to convert Rydberg to Hartree, and stores the
pure reference value in result metadata. It does not use a displaced energy or
add ZPE. For calculations using electronic smearing, the QE total energy can
contain a smearing contribution; isolated-molecule inputs should avoid
artificial electronic-temperature contributions.

A failed or interrupted calculation remains available for diagnosis. The next
calculation invocation creates a fresh numbered directory and seeds it from the immutable
reference density again.

Generated `scf.in` and `pw.in` files are self-contained for diagnosis: FDVIB
starts QE inside each attempt directory, uses `outdir='./out'`, rewrites
relative `pseudo_dir` paths for that location, and redirects `wfcdir` to the
same isolated `./out`. Users can enter an attempt directory and run
`pw.x -inp pw.in` directly. Treat completed attempt directories as immutable
snapshots; copy one before manual experiments. On the next invocation, FDVIB
may recover a valid manually completed attempt that has not yet been committed
to a completion snapshot.

Because `pw_command` is executed from the calculation directory, executables and
auxiliary files named in that command must either be discoverable through the
environment (for example `pw.x` on `PATH`) or use absolute paths.

FDVIB validates the total-force block, ignoring later decomposed force blocks
printed by some QE features. A successful task requires exactly `nat` total
forces and `JOB DONE` without an SCF non-convergence diagnostic.

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
and rotations are handled later by the appropriate analysis model.

With `run_dynmat=true`, FDVIB runs `dynmat_command` in an isolated attempt,
parses all `3N` frequencies and eigenvectors, and publishes:

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
`metadata.dat` is missing, exporters use conservative defaults:
`program=qe`, `mode_selection=all`, `multiplicity=1`, and
`electronic_energy_hartree=0.0`.

## State and recovery

`fdvib/state/dataset.state` fingerprints the SCF input and scientific dataset.
Changing geometry, selected atoms, displacement, system type, multiplicity, or
result prefix causes a hard error; use a different `outdir` for a new dataset.
Execution commands and `run_dynmat` are not part of that immutable identity.

The calculation holds an operating-system lock at `<outdir>.lock`, preventing
two FDVIB processes from writing the same dataset concurrently. The lock is
released automatically when the process exits.

Every completed stage has a commit marker and recorded file digests. The
restart behavior is:

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

The command reads `<prefix>.dynG`, `<prefix>.freq.out`, and optional
`metadata.dat`, then writes a CP2K-style, MfakeG-compatible `<prefix>.mol`.
It preserves imaginary-frequency signs, includes the complete geometry, and
writes zero displacement for frozen atoms. Re-running the command overwrites
the generated `.mol` file.

For `mode_selection=gas`, FDVIB writes only the molecular internal modes and
omits the artificial simulation cell. For `local` and `all`, it writes every
frequency exactly unequal to `0.0` and includes the three cell vectors in the
Molden `[Cell]` section. Atom and `FR-COORD` coordinates are written in Bohr,
while `[Cell]` vectors are written in Angstrom, matching the CP2K layout used
by MfakeG. The `[INT]` section is omitted because FDVIB does not calculate IR
intensities.

QE can also write Molden-style vibration files, but imaginary frequencies may
be inconvenient to inspect because some viewers treat them as zero-frequency
modes. FDVIB's Molden export keeps the signed frequency values.

### Difference from `dynmat.x` `filmol`

QE's `dynmat.x` has a `filmol` output whose default filename is
`dynmat.mold`; see the
[QE `INPUT_DYNMAT` reference](https://gitlab.com/QEF/q-e/blob/master/PHonon/Doc/INPUT_DYNMAT.txt).
FDVIB's `<prefix>.mol` is deliberately not a copy of that file:

| Property | `dynmat.x` `filmol` | FDVIB `<prefix>.mol` |
|---|---|---|
| Normal-mode data | Written directly by `dynmat.x` | Read from the same `dynmat.x` `filout` result |
| Frequency sign | Imaginary modes may appear as `0.00` in the `filmol` output | Signed `cm^-1` values from `filout` are preserved |
| Sections | Basic Molden frequency, coordinate, and normal-mode sections | Adds CP2K-style `[Cell]` and `[Atoms] AU` sections |
| Cell | No MfakeG cell metadata | Written for `local` and `all`; omitted for `gas` |
| Mode selection | QE writes the modes produced by `dynmat.x` | Gas internal modes, or all exactly nonzero modes for `local` and `all` |
| IR intensity | Not part of the basic `filmol` mode data | `[INT]` is omitted rather than fabricated |
| Main purpose | General Molden visualization | MfakeG-compatible conversion for GaussView |

FDVIB does not diagonalize a second dynamical matrix for this export. Both
files originate from the same `dynmat.x` calculation and use its normalized
phonon displacements. FDVIB reads `filout`, applies the metadata-controlled
mode selection, restores the cell and atomic metadata from `.dynG`, and writes
the CP2K-style Molden layout.

## GaussView visualization through MfakeG

After generating the Molden file, convert it with
[MfakeG](http://sobereva.com/soft/MfakeG):

```bash
fdvib modes fdvib/results
MfakeG fdvib/results/system.mol
```

MfakeG produces the Gaussian-like visualization file that GaussView can load
for frequency inspection and vibration animation. FDVIB does not generate or
maintain that Gaussian-like format. Because the `.mol` file has no `[INT]`
section, the converted visualization contains no calculated IR intensities.

## Shermo `.shm` export (recommended)

**Shermo is the recommended mature and comprehensive thermochemistry
post-processor for the frequency results.**
After `dynmat.x` has produced a valid frequency file, run:

```bash
fdvib shm fdvib/results
```

The `.shm` export produces an input file for
[Shermo](http://sobereva.com/soft/shermo).
> *Quoted from http://sobereva.com/soft/shermo:* Shermo is a free, general,
> very easy-to-use and flexible code for calculating molecular thermochemistry
> data based on ideal gas assumption.

> **If Shermo is utilized in your work, the following paper must be cited:**
>
> Tian Lu, Qinxue Chen, *Shermo: A general code for calculating molecular
> thermodynamic properties*, Comput. Theor. Chem., 1200, 113249 (2021)
> DOI: [10.1016/j.comptc.2021.113249](https://doi.org/10.1016/j.comptc.2021.113249)

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

For gas-molecule `.shm` files, run Shermo directly according to your target
analysis and the Shermo documentation, for example:

```bash
Shermo fdvib/results/molecule.shm
```

Export is transactional: FDVIB validates the generated text before renaming
it to `.shm`. Re-running the command overwrites the generated `.shm` file.

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
`mode_selection = gas` so the spin multiplicity is unambiguous. Results are written to
`fdvib/results/thermo.dat`. The file reports each energy table in eV,
kcal/mol, and kJ/mol; the corresponding entropy units are eV/K,
kcal/(mol K), and kJ/(mol K).

Changing temperature or other thermochemistry settings requires only another
`fdvib thermo` command; it never reruns the electronic-structure calculation.

## Acknowledgements

FDVIB uses Quantum ESPRESSO as its electronic-structure backend. Its `.shm`
interface follows the format and molecular thermochemistry workflow of
[Shermo](http://sobereva.com/soft/shermo), developed by Dr. Tian Lu at the
Beijing Kein Research Center for Natural Sciences and described in the work by
Tian Lu and Qinxue Chen cited above. We thank Dr. Tian Lu for answering a
question about representing molecular vibrational frequencies in Shermo
`.shm` files.

FDVIB's MfakeG-compatible Molden export follows the visualization workflow
demonstrated by Dr. Tian Lu's [MfakeG](http://sobereva.com/soft/MfakeG).
FDVIB independently writes the `.mol` file directly from FDVIB/QE results;
MfakeG can then convert it for GaussView.
