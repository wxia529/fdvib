#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Wanting Xia

set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'rm -rf "$case_dir"' EXIT

cat > "$case_dir/molecule.dynG" <<'EOF'
Dynamical matrix file
test molecule
  2    3   0  10.0000000 0 0 0 0 0
Basis vectors
  1 0 0
  0 1 0
  0 0 1
  1 'O      ' 14583.107889672
  2 'H      '  918.736996409
    1    1 0.0000000000 0.0000000000 0.0000000000
    2    2 0.1800000000 0.0000000000 0.0000000000
    3    2 0.0000000000 0.1700000000 0.0000000000
EOF

cat > "$case_dir/metadata.dat" <<'EOF'
program = qe
mode_selection = gas
multiplicity = 1
electronic_energy_hartree = -76.0
EOF

cat > "$case_dir/dynmat.in" <<'EOF'
asr = 'no'
remove_interaction_blocks = .false.
EOF

cat > "$case_dir/thermo.in" <<'EOF'
model = 'gas_rrho'
temperature_k = 298.15
pressure_atm = 1.0
symmetry_number = 1
rotor_type = 'nonlinear'
low_frequency_model = 'harmonic'
EOF

append_mode() {
  local index="$1" frequency="$2"
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$frequency" >> "$case_dir/molecule.freq.out"
  printf ' ( 1 0 0 0 0 0 )\n ( 0 0 1 0 0 0 )\n ( 0 0 0 0 1 0 )\n' >> "$case_dir/molecule.freq.out"
}

: > "$case_dir/molecule.freq.out"
append_mode 1 -0.3
append_mode 2 -0.2
append_mode 3 -0.1
append_mode 4 0.1
append_mode 5 0.2
append_mode 6 0.3
append_mode 7 100.0
append_mode 8 200.0
append_mode 9 300.0

"$fdvib" shm "$case_dir" > "$case_dir/shm.out"
grep -q 'SHM mode selection: nonlinear, retained 3, removed 6' "$case_dir/shm.out"
awk '/^\*wavenum/{inside=1;next}/^\*atoms/{inside=0}inside{print}' "$case_dir/molecule.shm" > "$case_dir/shm-frequencies"
test "$(wc -l < "$case_dir/shm-frequencies")" -eq 3
grep -qx '100.0000000000' "$case_dir/shm-frequencies"
grep -qx '200.0000000000' "$case_dir/shm-frequencies"
grep -qx '300.0000000000' "$case_dir/shm-frequencies"

"$fdvib" thermo "$case_dir" -inp "$case_dir/thermo.in" > /dev/null
grep -q '^# rotor_type: nonlinear$' "$case_dir/thermo.dat"
grep -q '^# rigid_body_modes_excluded: 6$' "$case_dir/thermo.dat"
grep -q '^# max_rigid_body_frequency_cm1: 0.3$' "$case_dir/thermo.dat"
grep -q '^# expected_vibrational_modes: 3$' "$case_dir/thermo.dat"
grep -q '^# positive_modes_used: 3$' "$case_dir/thermo.dat"
grep -q '^# units: T=K energies=kcal/mol entropy=kcal/mol/K$' "$case_dir/thermo.dat"
grep -q '^# units: T=K energies=kJ/mol entropy=kJ/mol/K$' "$case_dir/thermo.dat"
awk '
  /^# units: T=K energies=eV / { getline; getline; ev=$2 }
  /^# units: T=K energies=kcal\/mol / { getline; getline; kcal=$2 }
  /^# units: T=K energies=kJ\/mol / { getline; getline; kj=$2 }
  END {
    if (!(kcal/ev > 23.06054 && kcal/ev < 23.06056 &&
          kj/ev > 96.48532 && kj/ev < 96.48534)) exit 1
  }
' "$case_dir/thermo.dat"

sed -i "s/rotor_type = 'nonlinear'/rotor_type = 'linear'/" "$case_dir/thermo.in"
"$fdvib" thermo "$case_dir" -inp "$case_dir/thermo.in" > /dev/null
grep -q '^# rigid_body_modes_excluded: 5$' "$case_dir/thermo.dat"
grep -q '^# expected_vibrational_modes: 4$' "$case_dir/thermo.dat"
grep -q '^# positive_modes_used: 4$' "$case_dir/thermo.dat"

sed -i "s/rotor_type = 'linear'/rotor_type = 'nonlinear'/" "$case_dir/thermo.in"
sed -i 's/= 300.0 \[cm-1\]/= -300.0 [cm-1]/' "$case_dir/molecule.freq.out"
if "$fdvib" thermo "$case_dir" -inp "$case_dir/thermo.in" > /dev/null 2> "$case_dir/error"; then
  echo "gas RRHO unexpectedly accepted an imaginary vibrational mode" >&2
  exit 1
fi
grep -q 'imaginary/non-positive vibrational mode' "$case_dir/error"

cat > "$case_dir/thermo.in" <<'EOF'
model = 'local_harmonic'
temperature_k = 298.15
low_frequency_model = 'frequency_floor'
zero_tolerance_cm1 = 1.0
EOF
sed -i 's/remove_interaction_blocks = .false./remove_interaction_blocks = .true./' "$case_dir/dynmat.in"
"$fdvib" thermo "$case_dir" -inp "$case_dir/thermo.in" > /dev/null
grep -q '^# imaginary_modes_excluded: 1$' "$case_dir/thermo.dat"
grep -q '^# zero_modes_excluded: 6$' "$case_dir/thermo.dat"
grep -q '^# positive_modes_used: 2$' "$case_dir/thermo.dat"
grep -q '^# frequency_floor_cm1: 100$' "$case_dir/thermo.dat"
