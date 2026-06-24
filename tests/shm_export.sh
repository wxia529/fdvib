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

cat > "$case_dir/fdvib.in.reference" <<'EOF'
system_type = 'gas'
selected_atoms = 'all'
multiplicity = 1
EOF

cat > "$case_dir/electronic_structure.dat" <<'EOF'
electronic_energy_hartree = -76.000000000000000
multiplicity = 1
source = 'scf.out'
EOF
cat > "$case_dir/dynmat.in" <<'EOF'
asr = 'no'
remove_interaction_blocks = .false.
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
append_mode 5 0.42
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

if "$fdvib" shm "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "SHM export unexpectedly overwrote an existing file" >&2
  exit 1
fi
grep -q 'Refuse to overwrite' "$case_dir/error"

sed -i "s/asr = 'no'/asr = 'zero-dim'/" "$case_dir/dynmat.in"
if "$fdvib" shm "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "SHM export unexpectedly accepted asr='zero-dim'" >&2
  exit 1
fi
grep -q "requires asr='no'" "$case_dir/error"

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
cat > "$atom_dir/fdvib.in.reference" <<'EOF'
system_type = 'gas'
selected_atoms = 'all'
multiplicity = 1
EOF
cat > "$atom_dir/electronic_structure.dat" <<'EOF'
electronic_energy_hartree = -2.9
multiplicity = 1
source = 'scf.out'
EOF
cat > "$atom_dir/dynmat.in" <<'EOF'
asr = 'no'
remove_interaction_blocks = .false.
EOF
: > "$atom_dir/he.freq.out"
for index in 1 2 3; do
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$index" >> "$atom_dir/he.freq.out"
  printf ' ( 1 0 0 0 0 0 )\n' >> "$atom_dir/he.freq.out"
done
"$fdvib" shm "$atom_dir" > "$atom_dir/shm.out"
grep -q 'SHM mode selection: atom, retained 0, removed 3' "$atom_dir/shm.out"
test "$(awk '/^\*wavenum/{inside=1;next}/^\*atoms/{inside=0}inside{count++}END{print count+0}' "$atom_dir/he.shm")" -eq 0
