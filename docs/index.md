# FDVIB reference

FDVIB is a finite-difference vibration controller for Quantum ESPRESSO. A
calculation runs from one command and can be resumed safely:

```bash
fdvib -in fdvib.in
```

Normal-mode visualization and thermochemistry remain separate analyses:

```bash
fdvib modes fdvib/results
fdvib thermo fdvib/results -in thermo.in
fdvib shm fdvib/results
```

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

A new calculation starts with:

```text
calculation/
|-- scf.in
|-- fdvib.in
`-- thermo.in       optional, used only by the separate thermo command
```

`scf.in` stays external; it is not embedded in `fdvib.in`.

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

FDVIB currently requires `ibrav=0`, Cartesian positions in Angstrom, and cell
vectors in Angstrom. `tprnfor=.true.` is required for force extraction.

Do not set this in `scf.in`:

```fortran
startingpot = 'file'
```

FDVIB rejects that value before starting a calculation. The unperturbed
reference SCF uses the starting-potential policy from `scf.in` (normally QE's
default atomic superposition). FDVIB injects `startingpot='file'` only into
the generated displaced inputs.

Every displaced attempt has an independent QE `outdir`. Before execution,
FDVIB copies the reference `charge-density.dat` or `charge-density.hdf5` into
that attempt's `<qe-prefix>.save` directory. Tasks never share a writable QE
scratch directory.

The original `scf.in` should be tested independently. Tight electronic
convergence is important because force noise directly affects frequencies.

## fdvib.in

Local periodic example:

```fortran
&FDVIB
  scf_input             = 'scf.in',
  outdir                = 'fdvib',
  system_type           = 'local',
  selected_atoms        = 65, 66, 67,
  displacement_angstrom = 0.01,
  pw_command            = 'mpirun -np 8 pw.x',
  prefix                = 'system',
  run_dynmat            = .true.,
  dynmat_command        = 'dynmat.x',
/
```

Gas molecule example:

```fortran
&FDVIB
  scf_input             = 'scf.in',
  outdir                = 'fdvib',
  system_type           = 'gas',
  selected_atoms        = 'all',
  displacement_angstrom = 0.01,
  multiplicity          = 1,
  pw_command            = 'mpirun -np 8 pw.x',
  prefix                = 'molecule',
  run_dynmat            = .true.,
  dynmat_command        = 'dynmat.x',
/
```

The parser accepts one `key=value` assignment per line, ignores `!` comments,
and rejects unknown parameters.

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

## Reference SCF and displacement attempts

The reference stage is stored under:

```text
fdvib/reference/
|-- attempt-001/
|   |-- scf.in
|   |-- scf.out
|   `-- out/<qe-prefix>.save/charge-density.*
`-- complete.state
```

A reference attempt is committed only after FDVIB verifies the process exit
status, `JOB DONE`, SCF convergence, the complete force block, and a nonempty
charge-density file. FDVIB reads the last `! total energy = ... Ry` record from
the unperturbed reference SCF, divides it by two to convert Rydberg to Hartree,
and stores the pure reference value in `electronic_structure.dat`. It does not
use a displaced energy or add ZPE. For calculations using electronic smearing,
the QE total energy can contain a smearing contribution; isolated-molecule
inputs should avoid artificial electronic-temperature contributions.

Each displacement uses append-only attempts:

```text
fdvib/jobs/disp_0065_x_p/
|-- attempt-001/
|   |-- pw.in
|   |-- pw.out
|   `-- out/<qe-prefix>.save/charge-density.*
|-- forces.dat
`-- complete.state
```

A failed or interrupted attempt remains available for diagnosis. The next
calculation invocation creates a fresh attempt and seeds it from the immutable
reference density again.

FDVIB validates the total-force block, ignoring later decomposed force blocks
printed by some QE features. A successful task requires exactly `nat` total
forces and `JOB DONE` without an SCF non-convergence diagnostic.

## Hessian and QE dynamical matrix

For the active atom set, FDVIB evaluates

$$
H_{j\beta,i\alpha} =
-\frac{F_{j\beta}(+\delta_{i\alpha})-
F_{j\beta}(-\delta_{i\alpha})}{2\delta}.
$$

It reports the maximum antisymmetric component before symmetrization, then
writes:

```text
fdvib/results/<prefix>.dynG
fdvib/results/dynmat.in
fdvib/results/fdvib.in.reference
fdvib/results/electronic_structure.dat
```

Local calculations set `remove_interaction_blocks=.true.`. Gas calculations
retain the complete molecular Hessian. Both use `asr='no'`; rigid translations
and rotations are handled later by the appropriate analysis model.

With `run_dynmat=.true.`, FDVIB runs `dynmat_command` in an isolated attempt,
parses all `3N` frequencies and eigenvectors, and publishes:

```text
fdvib/results/dynmat.out
fdvib/results/<prefix>.freq.out
```

With `run_dynmat=.false.`, calculation stops after `dynG` and `dynmat.in`.
Changing only this setting to `.true.` and repeating `fdvib -in fdvib.in`
skips reference, displacement, and Hessian stages and runs only `dynmat.x`.

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

The command reads `<prefix>.dynG` and `<prefix>.freq.out`, then writes a compact
`<prefix>.mold`. It preserves imaginary-frequency signs, includes the complete
geometry, writes zero displacement for frozen atoms, and omits modes within
the internal `1e-6 cm^-1` zero tolerance.

## Thermochemistry

Thermochemistry is deliberately separate from the expensive calculation:

```bash
fdvib thermo fdvib/results -in thermo.in
```

Local harmonic example:

```fortran
&THERMO
  model                  = 'local_harmonic',
  temperature_k          = 298.15,
  low_frequency_model    = 'frequency_floor',
  frequency_floor_cm1    = 50.0,
  zero_tolerance_cm1     = 1.0,
/
```

Gas RRHO example:

```fortran
&THERMO
  model                    = 'gas_rrho',
  temperature_k            = 298.15,
  pressure_atm             = 1.0,
  symmetry_number          = 1,
  electronic_degeneracy    = 'auto',
  rotor_type               = 'auto',
  low_frequency_model      = 'harmonic',
/
```

Local analysis reports vibrational ZPE, internal energy, entropy, and free
energy. Gas RRHO additionally includes translation, rotation, and electronic
degeneracy. It removes three rigid modes for an atom, five for a linear
molecule, or six for a nonlinear molecule, then rejects any remaining
non-positive vibrational frequency. Results are written to
`fdvib/results/thermo.dat`.

Changing temperature or other thermochemistry settings requires only another
`fdvib thermo` command; it never reruns the electronic-structure calculation.

## Shermo `.shm` export

After `dynmat.x` has produced a valid frequency file, run:

```bash
fdvib shm fdvib/results
```

The command writes `<prefix>.shm` with exactly the four sections expected by
Shermo 2.6.2: electronic energy, wavenumbers, atoms, and electronic levels.
The electronic energy is the unperturbed reference SCF energy in Hartree;
coordinates are written in Angstrom, masses in amu, and the ground-state
electronic degeneracy is the configured multiplicity.
See the [Shermo 2.6.2 SHM compatibility specification](../SHM_COMPATIBILITY_SPEC_ZH.md)
for the exact text grammar and parser constraints.
QE species labels such as `C1` or `O_ads` are reduced to their leading
standard element symbol; labels that cannot be mapped to H--Og are rejected.
The exporter requires `asr='no'` and cross-checks
`remove_interaction_blocks` against the saved gas/local dataset before using
the frequencies.

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

For isolated gas molecules, use Shermo's normal molecular mode:

```bash
Shermo fdvib/results/molecule.shm
```

Export is transactional: FDVIB validates the generated text before renaming
it to `.shm` and refuses to overwrite an existing file.
