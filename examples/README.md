# FDVIB examples

The examples use `scf.in` as the QE input filename, but `scf_input` may point
to any filename. First make sure that `pw.x` can run your QE input successfully.
Then copy the appropriate `fdvib.in` and `thermo.in` into the same directory
and set the atom indices, execution commands, and other options required for
your system. The QE input must not set `startingpot='file'`; FDVIB runs an
unperturbed reference SCF and manages the seeded charge densities itself.

- `local/`: local vibrations for selected atoms in a periodic environment. The
  thermochemistry example uses a $100\ \mathrm{cm}^{-1}$ frequency floor.
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
