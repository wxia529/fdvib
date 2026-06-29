# Theory and method

This page describes the matrices and selection rules used by FDVIB. The main
reference page focuses on how to run calculations.

## Overview

If `dynmat.x` is run, FDVIB obtains vibrational modes in three conceptual
steps. First, finite differences of QE forces produce the Cartesian
force-constant matrix $\mathbf{C}$, which is the Cartesian Hessian of the
energy. Second, FDVIB writes that matrix in a QE-compatible Gamma-point
dynamical-matrix file. Third, `dynmat.x` performs mass weighting,
diagonalization, and frequency conversion.

FDVIB then selects the resulting modes for visualization, SHM export, or
thermochemistry. With `run_dynmat=false`, the workflow stops after writing
`dynG` and `dynmat.in`.

```text
QE forces -> force-constant matrix C -> dynG -> dynmat.x -> modes -> analysis
```

FDVIB does not diagonalize the force-constant matrix independently.

## Coordinate indexing and units

For $N$ atoms, FDVIB uses a $3N$ Cartesian coordinate vector. In the notation
below, $i$ and $j$ are one-based atom indices. The direction indices $\alpha$
and $\beta$ take the values 0, 1, and 2 for $x$, $y$, and $z$, respectively.
The corresponding combined indices are

$$
\begin{aligned}
p &= 3(i-1)+\alpha, \\
q &= 3(j-1)+\beta.
\end{aligned}
$$

Atom indices therefore match the one-based numbering used by QE and
`selected_atoms`. Direction indices are internal labels used only in the
matrix notation.

The displacement amplitude configured by `displacement_angstrom` is converted
to Bohr before force derivatives are formed:

$$
\delta_{\mathrm{Bohr}}
= \frac{\mathtt{displacement\_angstrom}}{\mathtt{BOHR\_TO\_ANG}}.
$$

QE reports forces in `Ry/Bohr`, so the displacement in the finite-difference
denominator must also be in Bohr. The resulting Cartesian force constants are
stored in `Ry/Bohr^2`.

## Central finite differences

For each selected atom $i$ and direction $\alpha$, FDVIB runs two SCF force
calculations:

$$
R_{i\alpha}^{(+)}=R_{i\alpha}+\delta,
\qquad
R_{i\alpha}^{(-)}=R_{i\alpha}-\delta.
$$

QE prints forces for all atoms. In a local calculation, FDVIB uses only the
force components on the selected atoms; forces on unselected environment atoms
are not included in the local force-constant matrix. The Cartesian
force-constant element is

$$
C_{j\beta,i\alpha}
=-\frac{F_{j\beta}(+\delta_{i\alpha})
-F_{j\beta}(-\delta_{i\alpha})}{2\delta}.
$$

The minus sign appears because the force-constant matrix is the second
derivative of the energy, while QE reports the negative energy gradient:

$$
\mathbf{F}=-\frac{\partial E}{\partial \mathbf{R}},
\qquad
\mathbf{C}=\frac{\partial^2 E}{\partial \mathbf{R}\,\partial \mathbf{R}}
=-\frac{\partial \mathbf{F}}{\partial \mathbf{R}}.
$$

## Force-constant matrix and symmetrization

FDVIB stores this as a dense $3N \times 3N$ matrix because a QE dynamical-matrix
file contains the full atomic structure:

$$
\begin{aligned}
\mathrm{row} &= 3(j-1)+\beta=q, \\
\mathrm{col} &= 3(i-1)+\alpha=p,
\end{aligned}
\qquad C[\mathrm{row},\mathrm{col}]=C_{j\beta,i\alpha}.
$$

For `system_type = gas`, `selected_atoms = all`, so every row and column is
computed. For `system_type = local`, only the selected-atom subblock contains
finite-difference force constants. All other blocks are zero placeholders used
to retain the full geometry in the QE-compatible file; they are not physical
force constants for the environment.

After all finite differences are read, FDVIB symmetrizes the matrix:

$$
C_{\mathrm{sym},pq}=\frac{C_{pq}+C_{qp}}{2}.
$$

For a local calculation, this operation effectively symmetrizes only the
selected-atom subblock because the other blocks are zero.

Before replacing the matrix with $C_{\mathrm{sym}}$, FDVIB reports the maximum
absolute antisymmetric component. A large value can indicate noisy forces,
insufficient SCF convergence, an unsuitable displacement, or an inconsistent
selected region.

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
rotations are handled later during analysis.

For local calculations, FDVIB sets:

```text
remove_interaction_blocks = .true.
```

For gas calculations, it sets:

```text
remove_interaction_blocks = .false.
```

This setting is used so that the zero-padded local matrix is handled
consistently during QE post-processing. It is disabled for gas calculations
because the complete molecular matrix is physical.

## Mass weighting and frequencies

The physical normal-mode problem is the mass-weighted eigenvalue equation:

$$
\sum_{j\beta}D_{i\alpha,j\beta}e_{k,j\beta}
=\lambda_k e_{k,i\alpha}.
$$

Here the mass-weighted dynamical matrix is

$$
D_{i\alpha,j\beta}
=\frac{C_{i\alpha,j\beta}}{\sqrt{M_iM_j}}.
$$

FDVIB does not diagonalize this matrix directly. It writes the QE-compatible
`dynG` file and lets `dynmat.x` perform the mass weighting, diagonalization,
frequency conversion, and eigenvector output.

Therefore, frequencies reported or exported by FDVIB are parsed from
`dynmat.x`; they do not come from a separate FDVIB diagonalization.

When `run_dynmat = true`, FDVIB runs `dynmat.x` and writes:

```text
fdvib/results/dynmat.out
fdvib/results/<prefix>.freq.out
```

FDVIB then parses exactly $3N$ modes from `<prefix>.freq.out`. Each mode
contains:

- a frequency in cm$^{-1}$;
- a three-component displacement vector for every atom.

Negative frequencies are preserved as negative values. They represent
imaginary modes in the usual vibrational-analysis convention.

The Molden-style `.mol` writer does not perform another diagonalization.
It combines these same `dynmat.x` frequencies and normalized displacements
with geometry and cell information parsed from `.dynG`.

## Local calculations

For `system_type = local`, FDVIB computes finite differences only for
`selected_atoms`. This is intended for a local harmonic description of an
active region embedded in an unselected environment.

The generated `<prefix>.dynG` still contains the full atomic geometry.
Unselected environment atoms therefore remain visible in Molden visualization
files, but their force-constant rows and columns outside the selected local
block are zero.

The resulting modes are a local harmonic approximation for the selected
region, not the complete vibrational spectrum of the periodic system. The
zero-padded environment blocks must not be interpreted as physical force
constants.

Post-processing then decides how many frequencies to keep:

| Export | Modes retained for `mode_selection = local` |
|---|---|
| SHM | $3N_{\mathrm{active}}$ modes |
| Molden visualization | All nonzero parsed modes |

This difference is intentional: SHM export is intended for thermochemical
use, while Molden export is intended for inspecting all modes produced by
`dynmat.x`.

For SHM local export, FDVIB removes the modes closest to zero until
$3N_{\mathrm{active}}$ modes remain, giving one mode for each Cartesian degree
of freedom in the active region. For example, selecting three atoms produces
nine modes in the SHM export.

## Gas calculations

For `system_type = gas`, FDVIB requires:

```text
selected_atoms = all
```

The full molecular force-constant matrix is constructed. The resulting $3N$
modes contain rigid translations, rigid rotations, and internal vibrations.
FDVIB removes rigid-body modes only during gas-mode post-processing.

The number of retained internal modes is

$$
N_{\mathrm{vib}}=
\begin{cases}
0, & \text{atom}, \\
3N-5, & \text{linear molecule}, \\
3N-6, & \text{nonlinear molecule}.
\end{cases}
$$

FDVIB classifies the molecule from the principal moments of inertia. A molecule
with one near-zero principal moment is treated as linear. It then sorts the
modes by absolute frequency, removes the expected 3, 5, or 6 rigid-body modes,
and restores the retained modes to their original order. The linearity test
uses an internal numerical tolerance.

The retained frequencies are not averaged. Numerically split degenerate modes
are preserved as separately calculated original frequencies. For example, the
two CO2 bending components are both written.

A true transition-state imaginary mode is retained when it is not among the
rigid modes closest to zero. If the largest removed absolute frequency and the
smallest retained absolute frequency are not clearly separated, inspect the
normal modes and improve geometry/force convergence before using the result for
thermochemistry.

## Export mode selection

`metadata.dat` controls how exporters select frequencies. These rules affect
post-processing only; they do not change the completed force calculations or
the force-constant matrix. See [Result metadata](index.md#result-metadata) for
the file format.

| `mode_selection` | SHM export | Molden visualization export |
|---|---|---|
| `gas` | Molecular internal modes | Molecular internal modes |
| `local` | $3N_{\mathrm{active}}$ modes | All nonzero modes |
| `all` | All nonzero modes | All nonzero modes |

Only frequencies stored as exactly `0.0` are omitted. Very small nonzero
frequencies are preserved.

## Thermochemistry

The built-in `thermo` command reads the parsed normal modes and evaluates the
model selected in `thermo.in`.

### Gas RRHO

With `model = gas_rrho`, `metadata.dat` must contain
`mode_selection = gas`. FDVIB removes rigid-body modes using the gas-mode rule
above and rejects imaginary or zero frequencies remaining in the vibrational
set. The rigid-rotor harmonic-oscillator model includes translational,
rotational, vibrational, and electronic-degeneracy contributions. The
electronic energy is the unperturbed SCF energy recorded in `metadata.dat`.
The rotational symmetry number is read from `thermo.in`; FDVIB does not infer
molecular point-group symmetry automatically. When
`electronic_degeneracy = auto`, the electronic degeneracy is taken from the
`multiplicity` field in `metadata.dat`.

### Local harmonic

With `model = local_harmonic`, FDVIB uses the selected local modes and applies
the configured low-frequency treatment. It includes only local vibrational
contributions; translation and rotation of the full periodic system are not
included.

`thermo.dat` reports thermochemical energies in eV, kcal/mol, and kJ/mol.
Entropies use the corresponding energy unit per kelvin.

## Infrared intensities

The finite-difference workflow described above constructs force constants and
normal modes. It does not compute:

- dipole derivatives;
- Born effective charges;
- dielectric tensors.

Therefore the Molden-style `.mol` export omits the `[INT]` section rather
than writing artificial zero intensities. Nonzero IR intensities require
additional electric-response information, which is outside the current FDVIB
finite-difference force workflow.
