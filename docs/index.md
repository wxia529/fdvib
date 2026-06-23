---
layout: default
title: FDVIB Reference
---

<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

# FDVIB Reference

FDVIB is a C++17 external finite-difference vibration controller for Quantum
ESPRESSO (QE). It creates positive and negative Cartesian displacements, runs
`pw.x` to obtain atomic forces, constructs a q=0 dynamical matrix, and uses
QE's `dynmat.x` for frequencies and normal modes.

FDVIB supports two calculation types:

- `local`: local vibrations in a periodic environment, such as an
  adsorbate or an active-site fragment. Atoms outside the selected region are
  treated as a frozen environment.
- `gas`: an isolated molecule in a vacuum cell, with ideal-gas rigid-rotor
  harmonic-oscillator (RRHO) thermochemical corrections.

FDVIB does not calculate full Brillouin-zone phonons and does not replace
`ph.x`, `q2r.x`, or Phonopy.

## Processing pipeline

```text
scf.in + fdvib.in
        |
        v
fdvib prepare
        |  create one directory for every atom/direction/sign
        v
fdvib run
        |  run pw.x, preserve output, extract all atomic forces
        v
fdvib analyze
        |  central differences -> Hessian -> system.dynG + dynmat.in
        v
dynmat.x
        |  create system.freq.out
        +----------------+
        |                |
        v                v
fdvib modes        fdvib thermo
        |                |
        v                v
system.mold        thermo.dat
```

Changing temperature, pressure, symmetry number, or the low-frequency model
requires only another `fdvib thermo` run. It does not require new `pw.x`
calculations.

## Requirements

- A C++17 compiler, such as GCC 8 or newer.
- CMake 3.12 or newer.
- Quantum ESPRESSO `pw.x`.
- Quantum ESPRESSO `dynmat.x`.

`pw.x` and `dynmat.x` may come from the same QE installation or may be made
available through an HPC module environment.

## Build and installation

From the FDVIB source directory:

```bash
cmake -S . -B build
cmake --build build -j2
```

The executable is created at:

```text
build/fdvib
```

The executable may be invoked from the build directory or installed:

```bash
cmake --install build --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

Tagged releases include a prebuilt Linux x86_64 package targeting glibc 2.17
or newer. `libstdc++` and `libgcc` are linked statically, so the package is
compatible with CentOS/RHEL 7 and most newer glibc distributions. Quantum
ESPRESSO executables are not included.

Verify the installation:

```bash
fdvib --version
fdvib --help
```

The command interface is:

```text
fdvib prepare fdvib.in [--force]
fdvib run fdvib.in
fdvib analyze fdvib.in
fdvib modes RESULTS_DIR
fdvib thermo RESULTS_DIR
```

`fdvib.in` and `thermo.in` use a simplified QE namelist style. Put only one
`key=value` assignment on each line and use `!` for comments. Do not combine
multiple assignments on one line.

## QE input requirements

The project directory must contain a validated `scf.in` with the following
properties:

```fortran
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
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

In particular:

- only `ATOMIC_POSITIONS angstrom` is accepted;
- only `CELL_PARAMETERS angstrom` is accepted;
- only `ibrav=0` is accepted;
- `tprnfor=.true.` is required so that QE prints forces;
- the unmodified `scf.in` should be tested before using FDVIB;
- force convergence directly controls frequency accuracy;
- a tight electronic `conv_thr` is recommended;
- every displaced calculation must use identical electronic settings.

Each generated input differs from `scf.in` in exactly two respects:

1. one Cartesian coordinate is displaced;
2. `outdir` points to the job-specific scratch directory.

All other QE settings are preserved.

## Project structure

A calculation directory has the following initial structure:

```text
my-calculation/
|-- scf.in
|-- fdvib.in
`-- thermo.in
```

Do not create displacement directories manually. `fdvib prepare` creates
them.

Configuration templates are provided in `examples/local/` and `examples/gas/`.
The following commands copy the local-periodic templates:

```bash
cp examples/local/fdvib.in .
cp examples/local/thermo.in .
```

## Local periodic calculations

Local mode analysis is intended for adsorbates, active sites, and local
functional groups. Example `fdvib.in`:

```fortran
&FDVIB
  scf_input               = 'scf.in',
  workdir                 = 'fdvib',
  system_type             = 'local',
  selected_atoms          = 65, 66, 67,
  displacement_angstrom   = 0.01,
  pw_command              = 'pw.x',
  mpi_command             = 'mpirun -np 8',
  output_prefix           = 'system',
/
```

| Parameter | Meaning |
|---|---|
| `scf_input` | QE SCF template |
| `workdir` | displacement jobs and result directory |
| `system_type` | `local` or `gas` |
| `selected_atoms` | one-based active atom indices |
| `displacement_angstrom` | positive/negative displacement magnitude in Angstrom |
| `pw_command` | `pw.x` command or path |
| `mpi_command` | launcher and MPI arguments; use an empty string for serial QE |
| `output_prefix` | result filename prefix |

For the active atom set $$A$$, FDVIB constructs the active Hessian from

$$
H_{j\beta,i\alpha} =
-\frac{F_{j\beta}(+\delta_{i\alpha})-
F_{j\beta}(-\delta_{i\alpha})}{2\delta}.
$$

Atoms outside the selected region remain frozen. The result is a frozen-
environment local vibration calculation, not a complete periodic phonon
spectrum.

## Gas-phase calculations

All atoms must be active for a gas molecule:

```fortran
&FDVIB
  scf_input               = 'scf.in',
  workdir                 = 'fdvib',
  system_type             = 'gas',
  selected_atoms          = 'all',
  displacement_angstrom   = 0.01,
  multiplicity            = 2,
  pw_command              = 'pw.x',
  mpi_command             = 'mpirun -np 8',
  output_prefix           = 'molecule',
/
```

For an open-shell molecule, the spin settings in `scf.in` must match the
multiplicity. A doublet requires, for example:

```fortran
&SYSTEM
  nspin = 2
  tot_magnetization = 1
/
```

A triplet requires:

```fortran
&SYSTEM
  nspin = 2
  tot_magnetization = 2
/
```

The relation for ordinary collinear spin is

$$
N_\uparrow-N_\downarrow=2S=M-1,
$$

where $$M$$ is the multiplicity. FDVIB validates this relation but does not edit
the spin settings in `scf.in`.

QE represents a gas molecule with a periodic cell. Accuracy therefore depends
on convergence with respect to vacuum size, periodic-image interactions,
dipole effects, and charged-cell corrections.

## `prepare` command

Run:

```bash
fdvib prepare fdvib.in
```

For three selected atoms, FDVIB creates

$$
3\ \text{atoms}\times3\ \text{directions}\times2\ \text{signs}=18
\ \text{jobs}.
$$

Example layout:

```text
fdvib/
|-- scf.in.reference
|-- fdvib.in.reference
|-- jobs.list
|-- disp_0065_x_p/
|   `-- disp_0065_x_p.in
|-- disp_0065_x_m/
|   `-- disp_0065_x_m.in
|-- disp_0065_y_p/
|   `-- disp_0065_y_p.in
`-- ...
```

Each displacement directory initially contains exactly one input file.

`scf.in.reference` and `fdvib.in.reference` are immutable snapshots of the
dataset definition. Later analysis reads the snapshots, not project files
that may have been edited after preparation.

The `fdvib.in` path passed to `analyze` is used only to locate `workdir`.
Geometry, selected atoms, displacement magnitude, system type, and output
prefix are loaded from the snapshots in that work directory.

### Repeated prepare and overwrite rules

- An identical existing input is preserved.
- A different existing input produces an error by default.
- Existing `pw.x` output is never overwritten.
- `--force` may update inputs only before any displacement output exists.
- FDVIB never deletes files automatically.

Before any calculation has started, inputs may be regenerated with:

```bash
fdvib prepare fdvib.in --force
```

## `run` command

Run:

```bash
fdvib run fdvib.in
```

FDVIB launches every `pw.x` calculation from the project root and redirects
standard output automatically. For example:

```text
fdvib/disp_0065_x_p/disp_0065_x_p.out
```

QE scratch data are written to:

```text
fdvib/disp_0065_x_p/out/
```

A completed job directory contains:

```text
disp_0065_x_p/
|-- disp_0065_x_p.in
|-- disp_0065_x_p.out
|-- forces.dat
|-- status
`-- out/
```

`forces.dat` contains the total forces from the last complete QE force block in
`Ry/Bohr`. Verbose non-local, ionic, Hubbard, SCF-correction, dispersion, and
other decomposed force contributions printed afterward are ignored.

FDVIB checks:

- the process exit code;
- the presence of `JOB DONE`;
- `convergence NOT achieved` diagnostics;
- `Error in routine` diagnostics;
- the presence of a complete force block;
- agreement between the number of forces and `nat`.

If an output file already exists, FDVIB will not rerun that job. It validates
the existing output and refreshes the derived `forces.dat` from its total-force
block. This also repairs force tables produced by older parser versions.

## `analyze` command

After all displacement jobs finish, run:

```bash
fdvib analyze fdvib.in
```

This step:

1. reads the positive and negative `forces.dat` files;
2. constructs the Hessian by central differences;
3. reports the maximum Hessian asymmetry before symmetrization;
4. symmetrizes the active Hessian;
5. writes a QE text dynamical matrix;
6. writes `dynmat.in`;
7. copies `fdvib.in.reference` into the result directory;
8. copies a project-level `thermo.in` into the result directory only when the
   destination does not already exist.

Output:

```text
fdvib/results/
|-- system.dynG
|-- fdvib.in.reference
|-- dynmat.in
`-- thermo.in
```

`analyze` does not execute `dynmat.x`.

If `thermo.in` did not exist during analysis, copy it later:

```bash
cp thermo.in fdvib/results/thermo.in
```

## QE post-processing

Enter the result directory and run QE post-processing:

```bash
cd fdvib/results
dynmat.x -in dynmat.in > dynmat.out
```

For local modes, the generated input is similar to:

```fortran
&INPUT
  fildyn='system.dynG',
  filout='system.freq.out',
  asr='no',
  remove_interaction_blocks=.true.,
/
```

Gas calculations also use `asr='no'`, but retain the complete molecular
Hessian:

```fortran
asr='no'
remove_interaction_blocks=.false.
```

Verify the output:

```bash
grep -E "freq|JOB DONE" dynmat.out
```

The result directory then contains `system.freq.out` and `dynmat.out`.
If `> dynmat.out` is omitted, `system.freq.out` is still written, but the run
log appears only on the terminal.

## `modes` command

Return to the project root and run:

```bash
fdvib modes fdvib/results
```

The command reads `system.dynG` and `system.freq.out`, then writes
`system.mold`.

Unlike QE's default Molden output, the compact FDVIB file:

- removes exact zero modes generated by the frozen environment;
- preserves the sign of imaginary frequencies;
- contains the complete geometry;
- writes zero displacement for frozen atoms;
- is suitable for Molden and Jmol.

This step is optional. It is not required for thermochemistry.

FDVIB refuses to overwrite an existing `system.mold`. Move or remove the old
file explicitly before regenerating it.

## Local harmonic configuration

Example `thermo.in` for a local calculation:

```fortran
&THERMO
  model                  = 'local_harmonic',
  temperature_k          = 298.15,
  low_frequency_model    = 'frequency_floor',
  frequency_floor_cm1    = 50.0,
  zero_tolerance_cm1     = 1.0,
/
```

Local thermochemistry does not use pressure, molecular translation, molecular
rotation or electronic degeneracy. FDVIB reports an input error if gas-only
parameters are supplied with `local_harmonic`.

## Gas RRHO configuration

Example `thermo.in` for a gas molecule:

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

Gas calculations require explicit values for:

- `pressure_atm`, in atmospheres (`1 atm = 101325 Pa`);
- `symmetry_number`, the rotational symmetry number.

The spin multiplicity is defined once in `fdvib.in`. `analyze` copies the
immutable configuration snapshot into the result directory, and `thermo`
reads the multiplicity from that snapshot.

With `electronic_degeneracy='auto'`,

$$
g_\mathrm{elec}=\text{multiplicity},
\qquad
S_\mathrm{elec}=k_B\ln g_\mathrm{elec}.
$$

If additional orbital degeneracy is known, provide an explicit positive
integer:

```fortran
electronic_degeneracy = 4
```

`rotor_type` accepts:

- `auto`: determine the type from principal moments of inertia;
- `atom`: monatomic species;
- `linear`: linear molecule;
- `nonlinear`: nonlinear molecule.

`auto` is appropriate unless a nearly linear geometry makes automatic
classification unstable.

The rotor type also fixes the vibrational degrees of freedom. FDVIB removes
the rigid-body modes closest to zero and requires exactly:

- `0` vibrational modes for an atom;
- `3N-5` vibrational modes for a linear molecule;
- `3N-6` vibrational modes for a nonlinear molecule.

FDVIB reports the largest absolute frequency among the removed rigid-body
modes as `max_rigid_body_frequency_cm1`. This is a diagnostic rather than a
hard cutoff: a clear separation between the removed modes and the retained
vibrations is more informative than a universal absolute threshold. Any
non-positive frequency remaining in the vibrational set is a fatal error:
optimize the geometry and repeat the frequency calculation before using gas
RRHO thermochemistry.

Rigid-body identification is frequency-based: FDVIB removes the required
number of modes with the smallest absolute frequencies. It does not project
the eigenvectors onto an explicit translational/rotational subspace. If a true
very soft vibration lies closer to zero than a residual rigid-body mode, the
classification can be ambiguous. Inspect the normal modes and improve force
convergence when the low-frequency ordering is uncertain.

## `thermo` command

Run:

```bash
fdvib thermo fdvib/results
```

FDVIB requires exactly one `*.dynG`, exactly one `*.freq.out`, and one
`thermo.in` in the result directory. Gas RRHO calculations additionally
require `fdvib.in.reference`. The command writes `thermo.dat`.

If `dynmat.in` is present, FDVIB also cross-checks the physical model:

- both models require `asr='no'`;
- `local_harmonic` requires `remove_interaction_blocks=.true.`;
- `gas_rrho` requires `remove_interaction_blocks=.false.`.

This prevents a local-mode frequency file from being used accidentally as a
gas-molecule calculation.

Example local output:

```text
# FDVIB thermochemistry
# model: local_harmonic
# low_frequency_model: frequency_floor
# frequency_floor_cm1: 50
# imaginary_modes_excluded: 3
# zero_modes_excluded: 192
# positive_modes_used: 6
# modes_floored: 1
# units: T=K energies=eV entropy=eV/K
#         T/K          ZPE/eV        U_vib/eV         S_vib/eV_K       TS_vib/eV        F_vib/eV
      298.150    0.3114394382    0.3649497141     0.000438489323    0.1307355917    0.2342141224
```

Machine-readable output uses only:

- temperature: K;
- frequency: `cm^-1`;
- energy: eV;
- entropy: eV/K.

Units such as kJ/mol and kcal/mol are not mixed into the same data table.

Gas output additionally records `rotor_type`,
`rigid_body_modes_excluded`, `max_rigid_body_frequency_cm1`, and
`expected_vibrational_modes`, so the `3N-5` or `3N-6` selection can be audited
directly.

## Thermochemical quantities

### ZPE

The zero-point vibrational energy is

$$
E_\mathrm{ZPE}=\frac12\sum_i h\nu_i.
$$

### U_vib

The vibrational internal energy includes ZPE and thermal excitation:

$$
U_\mathrm{vib}(T)=\sum_i\left[
\frac12h\nu_i+\frac{h\nu_i}{e^{h\nu_i/k_BT}-1}
\right].
$$

### S_vib and TS_vib

`S_vib` is harmonic vibrational entropy and

$$
TS_\mathrm{vib}=T S_\mathrm{vib}.
$$

### F_vib

The local vibrational Helmholtz free energy is

$$
F_\mathrm{vib}=U_\mathrm{vib}-TS_\mathrm{vib}.
$$

A local correction is commonly combined as

$$
E_\mathrm{corrected}(T)=E_\mathrm{DFT}+F_\mathrm{vib}(T).
$$

FDVIB does not read or combine the QE electronic energy automatically.

### Gas H_corr and G_corr

Gas RRHO thermochemistry also includes translation, rotation, and electronic
entropy:

$$
H_\mathrm{corr}=H_\mathrm{trans}+U_\mathrm{rot}+U_\mathrm{vib},
$$

$$
G_\mathrm{corr}=H_\mathrm{corr}-T\left(
S_\mathrm{trans}+S_\mathrm{rot}+S_\mathrm{vib}+S_\mathrm{elec}
\right).
$$

Here `H_trans = 3/2 k_BT + PV = 5/2 k_BT` for one ideal-gas molecule, so this
is equivalently

$$
G_\mathrm{corr}=E_\mathrm{ZPE}+\Delta U(0\rightarrow T)+PV-TS.
$$

Add `G_corr` to the electronic energy calculated at the same theoretical
level to obtain the ideal-gas Gibbs energy at the requested temperature and
pressure. FDVIB does not combine those two values automatically.

## Low-frequency models

FDVIB supports two models.

### harmonic

```fortran
low_frequency_model = 'harmonic'
```

Every positive frequency is used without modification. This is the default in
the gas example and is mandatory for `gas_rrho`.

### frequency_floor

```fortran
low_frequency_model = 'frequency_floor'
frequency_floor_cm1 = 50.0
```

For a positive frequency below the threshold,

$$
\nu_\mathrm{used}=\max(\nu_\mathrm{raw},\nu_\mathrm{floor}).
$$

The substituted frequency is used consistently for ZPE, `U_vib`, `S_vib`,
and free energy. The local example uses a `50 cm^-1` floor. This model is
rejected for gas RRHO calculations; molecular translation and rotation are
removed using the rigid-body degree count instead.

## Imaginary and zero modes

For `local_harmonic`, the fixed rules are:

- if `abs(frequency) < zero_tolerance_cm1`, classify and exclude it as zero;
- if the frequency is negative, report and exclude it as imaginary;
- process positive frequencies with the selected low-frequency model.

FDVIB never replaces an imaginary frequency by its absolute value in the
partition function, because an imaginary mode is not a stable harmonic
oscillator.

`thermo.dat` records:

```text
imaginary_modes_excluded
zero_modes_excluded
positive_modes_used
modes_floored
```

Interpretation of an imaginary mode must distinguish a physical instability
from insufficient geometry optimization or finite-difference noise.

For `gas_rrho`, FDVIB first removes exactly three atomic translations, five
linear-molecule rigid motions, or six nonlinear-molecule rigid motions. The
gas model does not use `zero_tolerance_cm1`; it reports the largest removed
absolute frequency for diagnosis. A negative or zero frequency in the
remaining `3N-5` or `3N-6` vibrational modes terminates the calculation; it is
not silently omitted or replaced by its absolute value.

## Recalculation at another temperature

Edit `fdvib/results/thermo.in`:

```fortran
temperature_k = 500.0
```

Then run:

```bash
fdvib thermo fdvib/results
```

This operation does not generate displacements, run `pw.x`, reconstruct the
Hessian, or run `dynmat.x`. It recalculates only partition functions and
thermochemical quantities.

To preserve several temperatures, copy the old table before the next run:

```bash
cp fdvib/results/thermo.dat fdvib/results/thermo_298.15K.dat
```

## Recalculation dependencies

| Changed item | thermo | dynmat.x | analyze | pw.x |
|---|---:|---:|---:|---:|
| temperature | yes | no | no | no |
| gas pressure | yes | no | no | no |
| rotational symmetry number | yes | no | no | no |
| electronic degeneracy | yes | no | no | no |
| low-frequency model or threshold | yes | no | no | no |
| Molden file | no | no | no | no |
| ASR in `dynmat.in` | yes | yes | no | no |
| displacement magnitude | yes | yes | yes | yes |
| selected atoms | yes | yes | yes | yes |
| geometry | yes | yes | yes | yes |
| QE convergence settings | yes | yes | yes | yes |

Regenerating only the Molden file requires `fdvib modes`, but the existing
`system.mold` must first be moved or removed explicitly.

## Data retention and cleanup

FDVIB does not implement automatic cleanup. After successful analysis, the QE
scratch directories may be removed manually:

```text
fdvib/disp_*/out/
```

The following files are required to preserve reproducibility:

```text
scf.in.reference
fdvib.in.reference
jobs.list
disp_*/disp_*.in
disp_*/disp_*.out
disp_*/forces.dat
results/system.dynG
results/system.freq.out
results/dynmat.in
results/fdvib.in.reference
results/thermo.in
```

As long as the displacement inputs, QE outputs, and force tables are retained,
the dynamical matrix can be reconstructed.

## Diagnostics

### Cannot find scf.in

Verify the path in `fdvib.in`:

```fortran
scf_input = 'scf.in'
```

It is resolved relative to the directory containing `fdvib.in`.

### Require ATOMIC_POSITIONS angstrom

Convert the structure explicitly to:

```text
ATOMIC_POSITIONS angstrom
```

Version 0.1 does not convert `crystal`, `alat`, or `bohr` automatically.

### scf.in must contain tprnfor=.true.

Add this setting to `&CONTROL`:

```fortran
tprnfor = .true.
```

Then rerun `prepare --force` before any displacement calculation has started.

### Existing input differs

The existing generated input does not match the current configuration. Before
any output exists, regenerate it with:

```bash
fdvib prepare fdvib.in --force
```

After any output exists, FDVIB refuses to mix dataset definitions. Use a new
work directory.

### JOB DONE not found

The corresponding `pw.x` calculation did not finish normally. Inspect its
output, correct the problem, move the failed output explicitly, and rerun.
FDVIB never overwrites the failed output.

### Incomplete force block

Confirm that:

- `tprnfor=.true.` is set;
- the QE output is not truncated;
- the output contains exactly `nat` force lines.

### Result directory must contain exactly one .dynG and one .freq.out

Run `dynmat.x` in the result directory:

```bash
dynmat.x -in dynmat.in > dynmat.out
```

Also move duplicate or backup `*.dynG` and `*.freq.out` files out of that
directory.

### Refuse to overwrite system.mold

FDVIB protects existing visualization output. Move or remove it explicitly,
then rerun `fdvib modes`.

### A local calculation contains many zero modes

This is expected. The QE dynamical matrix retains the full `3N` dimension,
while rows and columns for frozen atoms are zero. `fdvib modes` removes exact
zero modes, and `fdvib thermo` excludes them from thermochemistry.

### Frequencies depend strongly on displacement magnitude

Compare, for example, 0.005, 0.01, and 0.02 Angstrom. Strong sensitivity may
indicate:

- insufficient SCF force convergence;
- a displacement too small relative to force noise;
- a displacement too large for the harmonic region;
- an insufficiently optimized reference geometry.

## Limitations

- No full periodic phonon dispersion or phonon density of states.
- No q-point mesh integration.
- No IR or Raman intensities.
- No quasi-RRHO or hindered-rotor model.
- No automatic cleanup of QE scratch data.
- Only `ibrav=0` and Angstrom coordinates are supported.
- QE electronic energies are not combined automatically with corrections.
- Rotational symmetry numbers are not inferred automatically.
- Imaginary frequencies are never converted to positive frequencies for
  thermochemistry.

## Command summary

```bash
# 1. Create displacement inputs
fdvib prepare fdvib.in

# 2. Run all pw.x jobs
fdvib run fdvib.in

# 3. Construct the QE dynamical matrix
fdvib analyze fdvib.in

# 4. Run QE post-processing
cd fdvib/results
dynmat.x -in dynmat.in > dynmat.out
cd ../..

# 5. Optional: create a compact Molden file
fdvib modes fdvib/results

# 6. Calculate thermochemical corrections
fdvib thermo fdvib/results
```
