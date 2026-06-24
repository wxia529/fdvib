# FDVIB examples

The examples are configuration templates. Copy the appropriate `fdvib.in`
and `thermo.in` into the directory containing a validated QE `scf.in`, then
set atom indices, executable commands, and molecular metadata. `scf.in` must
not set `startingpot='file'`; FDVIB runs an unperturbed reference SCF and
manages the seeded charge densities itself.

- `local/`: frozen-environment local vibrations. The thermochemistry example
  uses a 50 `cm^-1` frequency floor.
- `gas/`: isolated-molecule RRHO thermochemistry. The thermochemistry example
  uses unmodified harmonic frequencies. Pressure in atm and the rotational
  symmetry number are set in `thermo.in`; multiplicity is set once in
  `fdvib.in`.

Example:

```bash
cp examples/local/fdvib.in .
cp examples/local/thermo.in .
fdvib -inp fdvib.in
fdvib modes fdvib/results
fdvib thermo fdvib/results -inp thermo.in
fdvib shm fdvib/results
```
