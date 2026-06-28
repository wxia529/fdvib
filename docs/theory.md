# Theory and method

This page describes the matrices and selection rules used by FDVIB. The main
reference page focuses on how to run calculations.

## Coordinate indexing

For `N` atoms, FDVIB uses a `3N` Cartesian coordinate vector. Atom indices in
input files are one-based, but the matrix notation below uses abstract indices:

```text
i, j    one-based atom indices
α, β    zero-based Cartesian directions: x=0, y=1, z=2
p       combined coordinate index p = 3*(i - 1) + α
q       combined coordinate index q = 3*(j - 1) + β
```

The displacement amplitude configured by `displacement_angstrom` is converted
to Bohr before force derivatives are formed:

```text
δ_bohr = displacement_angstrom / BOHR_TO_ANG
```

QE forces are read in `Ry/Bohr`, so the finite-difference force constants are
stored in `Ry/Bohr^2`.

## Force-constant matrix

For each active atom `i` and direction `α`, FDVIB runs two SCF force
calculations:

```text
R(iα) + δ
R(iα) - δ
```

It reads the total force on each active atom `j` and direction `β` from the two
outputs. The Cartesian force-constant element is then:

```text
C[jβ,iα] = - (F[jβ](+δ[iα]) - F[jβ](-δ[iα])) / (2δ)
```

The minus sign appears because the Hessian/force-constant matrix is the
derivative of the energy, while QE reports forces:

```text
F = - dE/dR
C =  d²E/dR dR = - dF/dR
```

Implementation-wise, FDVIB stores this as a dense `3N × 3N` matrix:

```text
C[row, col]
row = 3*(j - 1) + β
col = 3*(i - 1) + α
```

Only active displacement columns are computed. For `system_type = gas`,
`selected_atoms = all`, so all `3N` columns are computed. For
`system_type = local`, only selected atoms are displaced, and only force rows
for selected atoms are filled. The rest of the dense matrix remains zero.

After all finite differences are read, FDVIB symmetrizes the matrix:

```text
C_sym[p,q] = (C[p,q] + C[q,p]) / 2
```

Before replacing the matrix with `C_sym`, FDVIB reports the maximum absolute
antisymmetric component. A large value is a practical warning for force noise,
insufficient SCF convergence, too large/small displacement, or an inconsistent
local active region.

## QE dynamical matrix file

FDVIB writes the symmetrized Cartesian force-constant matrix to:

```text
fdvib/results/<prefix>.dynG
```

The file uses a QE Gamma-point dynamical-matrix layout:

- cell vectors are written from the original QE input;
- atomic positions are written in the QE dynamical-matrix convention;
- species masses are converted to QE Rydberg atomic-unit masses;
- the Cartesian matrix block is written as real values with zero imaginary
  parts.

FDVIB also writes:

```text
fdvib/results/dynmat.in
```

with:

```text
asr = 'no'
```

FDVIB does not ask `dynmat.x` to remove acoustic modes. Rigid translations and
rotations are handled later by the relevant post-processing model.

For local calculations, FDVIB sets:

```text
remove_interaction_blocks = .true.
```

For gas calculations, it sets:

```text
remove_interaction_blocks = .false.
```

The gas matrix is already a full molecular matrix. The local matrix contains
the selected local block embedded in the full geometry.

## Mass weighting and frequencies

The physical normal-mode problem is the mass-weighted eigenvalue equation:

```text
Σ[jβ] D[iα,jβ] e_k[jβ] = λ_k e_k[iα]
```

where the mass-weighted dynamical matrix is conceptually:

```text
D[iα,jβ] = C[iα,jβ] / sqrt(M_i M_j)
```

FDVIB does not diagonalize this matrix directly. It writes the QE-compatible
`dynG` file and lets `dynmat.x` perform the mass weighting, diagonalization,
frequency conversion, and eigenvector output.

When `run_dynmat = true`, FDVIB runs `dynmat.x` and publishes:

```text
fdvib/results/dynmat.out
fdvib/results/<prefix>.freq.out
```

FDVIB then parses exactly `3N` modes from `<prefix>.freq.out`. A mode contains:

```text
frequency in cm^-1
3-component displacement vector for each atom
```

Negative frequencies are preserved as negative values. They represent
imaginary modes in the usual vibrational-analysis convention.

The MfakeG-compatible `.mol` writer does not perform another diagonalization.
It combines these same `dynmat.x` frequencies and normalized displacements
with geometry and cell information parsed from `.dynG`.

## Local calculations

For `system_type = local`, FDVIB computes finite differences only for
`selected_atoms`. This is intended for a local harmonic description of an
active region embedded in a frozen environment.

The generated `<prefix>.dynG` still contains the full atomic geometry. Frozen
atoms therefore remain visible in Molden/visualization exports, but their
force-constant rows and columns outside the selected local block are zero.

Post-processing then decides how many frequencies to keep:

```text
SHM export:
  mode_selection = local   keep 3 * N_active modes

Molden visualization export:
  mode_selection = local   keep all nonzero parsed frequencies
```

For SHM local export, FDVIB removes the modes closest to zero until
`3 * N_active` modes remain. This is a practical compact representation for
Shermo local-mode input.

## Gas calculations

For `system_type = gas`, FDVIB requires:

```text
selected_atoms = all
```

The full molecular Hessian is constructed. The resulting `3N` modes contain
rigid translations, rigid rotations, and internal vibrations. FDVIB removes
rigid-body modes only during gas-mode post-processing.

The number of retained internal modes is:

```text
atom       0
linear     3N - 5
nonlinear  3N - 6
```

FDVIB classifies the molecule from the principal moments of inertia. A molecule
with one near-zero principal moment is treated as linear. Rigid-body modes are
identified as the modes with the smallest absolute frequencies:

```text
sort modes by |frequency|
remove 3, 5, or 6 lowest-|frequency| modes
retain the remaining modes in original order
```

The retained frequencies are not averaged. Numerically split degenerate modes
are preserved as separately calculated original frequencies. For example, the
two CO2 bending components are both written.

A true transition-state imaginary mode is retained when it is not among the
rigid modes closest to zero. If the largest removed absolute frequency and the
smallest retained absolute frequency are not clearly separated, inspect the
normal modes and improve geometry/force convergence before using the result for
thermochemistry.

## Export mode selection

`metadata.dat` controls how exporters select frequencies.

For SHM export:

```text
mode_selection = gas     molecular internal modes
mode_selection = local   3 * selected_atoms modes
mode_selection = all     all nonzero modes
```

For Molden visualization export:

```text
mode_selection = gas     molecular internal modes
mode_selection = local   all nonzero modes
mode_selection = all     all nonzero modes
```

Here “nonzero” means the parsed frequency is exactly not `0.0`. Very small but
nonzero frequencies are preserved.

## Thermochemistry

The built-in `thermo` command reads the parsed normal modes and evaluates
quantities according to `thermo.in`.

For:

```text
model = gas_rrho
```

FDVIB requires `metadata.dat` with:

```text
mode_selection = gas
```

It removes rigid-body modes using the gas-mode rule above and rejects
imaginary or zero frequencies remaining in the true vibrational set. The
electronic energy is the unperturbed reference SCF energy recorded in
`metadata.dat`.

For:

```text
model = local_harmonic
```

FDVIB uses the selected local normal modes and applies the configured
low-frequency treatment.

`thermo.dat` reports the same thermochemical quantities in eV, kcal/mol, and
kJ/mol. Entropies use the corresponding energy unit per kelvin.

## Infrared intensities

The finite-difference workflow described above constructs force constants and
normal modes. It does not compute:

- dipole derivatives;
- Born effective charges;
- dielectric tensors.

Therefore the MfakeG-compatible `.mol` export omits the Molden `[INT]` section
instead of inventing zero intensities. Nonzero IR intensities require
additional electric-response information, which is outside the current FDVIB
finite-difference force workflow.
