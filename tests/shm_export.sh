#!/usr/bin/env bash
set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'rm -rf "$case_dir"' EXIT

cat > "$case_dir/co2.dynG" <<'EOF'
Dynamical matrix file
CO2 test
  2    3   0  10.0000000 0 0 0 0 0
Basis vectors
  1 0 0
  0 1 0
  0 0 1
  1 'C1     ' 10937.3309176688
  2 'O1     ' 14582.1922119821
    1    2 -0.2200000000 0.0000000000 0.0000000000
    2    1  0.0000000000 0.0000000000 0.0000000000
    3    2  0.2200000000 0.0000000000 0.0000000000
EOF

cat > "$case_dir/metadata.dat" <<'EOF'
program = qe
mode_selection = gas
electronic_energy_hartree = -76.000000000000000
multiplicity = 1
EOF

append_mode() {
  local index="$1" frequency="$2"
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$frequency" >> "$case_dir/co2.freq.out"
  printf ' ( 1 0 0 0 0 0 )\n ( 0 0 1 0 0 0 )\n ( 0 0 0 0 1 0 )\n' >> "$case_dir/co2.freq.out"
}

: > "$case_dir/co2.freq.out"
append_mode 1 -25.66
append_mode 2 -16.95
append_mode 3 -6.30
append_mode 4 -0.77
append_mode 5 0.0
append_mode 6 630.07
append_mode 7 641.78
append_mode 8 1326.44
append_mode 9 2383.22

"$fdvib" shm "$case_dir" > "$case_dir/shm.out"
grep -q 'SHM mode selection: linear, retained 4, removed 5' "$case_dir/shm.out"
awk '/^\*wavenum/{inside=1;next}/^\*atoms/{inside=0}inside{print}' "$case_dir/co2.shm" > "$case_dir/frequencies"
test "$(wc -l < "$case_dir/frequencies")" -eq 4
grep -qx '630.0700000000' "$case_dir/frequencies"
grep -qx '641.7800000000' "$case_dir/frequencies"
grep -qx '1326.4400000000' "$case_dir/frequencies"
grep -qx '2383.2200000000' "$case_dir/frequencies"
grep -q '^-7.600000000000000E+01$' "$case_dir/co2.shm"
grep -q '^0.000000 1$' "$case_dir/co2.shm"
test "$(grep -Ec '^(C|O) ' "$case_dir/co2.shm")" -eq 3

"$fdvib" shm "$case_dir" > "$case_dir/shm-rerun.out"
grep -q 'SHM mode selection: linear, retained 4, removed 5' "$case_dir/shm-rerun.out"

rm "$case_dir/co2.shm"
sed -i 's/program = qe/program = vasp/' "$case_dir/metadata.dat"
if "$fdvib" shm "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "SHM export unexpectedly accepted unsupported program" >&2
  exit 1
fi
grep -q 'Unsupported program in metadata.dat: vasp' "$case_dir/error"
sed -i 's/program = vasp/program = qe/' "$case_dir/metadata.dat"

sed -i 's/mode_selection = gas/mode_selection = unknown/' "$case_dir/metadata.dat"
if "$fdvib" shm "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "SHM export unexpectedly accepted invalid mode_selection" >&2
  exit 1
fi
grep -q 'mode_selection must be all, gas, or local in metadata.dat' "$case_dir/error"
sed -i 's/mode_selection = unknown/mode_selection = gas/' "$case_dir/metadata.dat"

rm "$case_dir/metadata.dat"
"$fdvib" shm "$case_dir" > "$case_dir/default-shm.out"
grep -q 'SHM mode selection: all, retained 8, removed 1' "$case_dir/default-shm.out"
grep -q '^0.000000000000000E+00$' "$case_dir/co2.shm"
awk '/^\*wavenum/{inside=1;next}/^\*atoms/{inside=0}inside{print}' "$case_dir/co2.shm" > "$case_dir/default-frequencies"
test "$(wc -l < "$case_dir/default-frequencies")" -eq 8
if grep -qx '0.0000000000' "$case_dir/default-frequencies"; then
  echo "SHM export unexpectedly kept an exact zero frequency in all mode" >&2
  exit 1
fi
rm "$case_dir/co2.shm"
cat > "$case_dir/metadata.dat" <<'EOF'
program = qe
mode_selection = local
selected_atoms = 1,1
EOF
if "$fdvib" shm "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "SHM export unexpectedly accepted duplicate selected_atoms" >&2
  exit 1
fi
grep -q 'Invalid/duplicate selected atom in metadata.dat' "$case_dir/error"

atom_dir="$case_dir/atom"
mkdir -p "$atom_dir"
cat > "$atom_dir/he.dynG" <<'EOF'
Dynamical matrix file
He test
  1    1   0  10.0000000 0 0 0 0 0
Basis vectors
  1 0 0
  0 1 0
  0 0 1
  1 'He     ' 3647.7423100000
    1    1 0.0000000000 0.0000000000 0.0000000000
EOF
cat > "$atom_dir/metadata.dat" <<'EOF'
program = qe
mode_selection = gas
electronic_energy_hartree = -2.9
multiplicity = 1
EOF
: > "$atom_dir/he.freq.out"
for index in 1 2 3; do
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$index" >> "$atom_dir/he.freq.out"
  printf ' ( 1 0 0 0 0 0 )\n' >> "$atom_dir/he.freq.out"
done
"$fdvib" shm "$atom_dir" > "$atom_dir/shm.out"
grep -q 'SHM mode selection: atom, retained 0, removed 3' "$atom_dir/shm.out"
test "$(awk '/^\*wavenum/{inside=1;next}/^\*atoms/{inside=0}inside{count++}END{print count+0}' "$atom_dir/he.shm")" -eq 0
