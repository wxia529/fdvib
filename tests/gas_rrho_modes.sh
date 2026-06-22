#!/usr/bin/env bash
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

cat > "$case_dir/fdvib.in.reference" <<'EOF'
system_type = 'gas'
multiplicity = 1
EOF

cat > "$case_dir/dynmat.in" <<'EOF'
asr = 'zero-dim'
remove_interaction_blocks = .false.
EOF

cat > "$case_dir/thermo.in" <<'EOF'
model = 'gas_rrho'
temperature_k = 298.15
pressure_bar = 1.0
symmetry_number = 1
rotor_type = 'nonlinear'
low_frequency_model = 'harmonic'
zero_tolerance_cm1 = 1.0
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

"$fdvib" thermo "$case_dir" > /dev/null
grep -q '^# rotor_type: nonlinear$' "$case_dir/thermo.dat"
grep -q '^# rigid_body_modes_excluded: 6$' "$case_dir/thermo.dat"
grep -q '^# expected_vibrational_modes: 3$' "$case_dir/thermo.dat"
grep -q '^# positive_modes_used: 3$' "$case_dir/thermo.dat"

sed -i "s/rotor_type = 'nonlinear'/rotor_type = 'linear'/" "$case_dir/thermo.in"
"$fdvib" thermo "$case_dir" > /dev/null
grep -q '^# rigid_body_modes_excluded: 5$' "$case_dir/thermo.dat"
grep -q '^# expected_vibrational_modes: 4$' "$case_dir/thermo.dat"
grep -q '^# positive_modes_used: 4$' "$case_dir/thermo.dat"

sed -i "s/rotor_type = 'linear'/rotor_type = 'nonlinear'/" "$case_dir/thermo.in"
sed -i 's/= 300.0 \[cm-1\]/= -300.0 [cm-1]/' "$case_dir/molecule.freq.out"
if "$fdvib" thermo "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "gas RRHO unexpectedly accepted an imaginary vibrational mode" >&2
  exit 1
fi
grep -q 'imaginary/non-positive vibrational mode' "$case_dir/error"

cat > "$case_dir/thermo.in" <<'EOF'
model = 'local_harmonic'
temperature_k = 298.15
low_frequency_model = 'harmonic'
zero_tolerance_cm1 = 1.0
EOF
sed -i "s/asr = 'zero-dim'/asr = 'no'/" "$case_dir/dynmat.in"
sed -i 's/remove_interaction_blocks = .false./remove_interaction_blocks = .true./' "$case_dir/dynmat.in"
"$fdvib" thermo "$case_dir" > /dev/null
grep -q '^# imaginary_modes_excluded: 1$' "$case_dir/thermo.dat"
grep -q '^# zero_modes_excluded: 6$' "$case_dir/thermo.dat"
grep -q '^# positive_modes_used: 2$' "$case_dir/thermo.dat"
