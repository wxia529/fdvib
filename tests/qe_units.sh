#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Wanting Xia

set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'if [[ -n "${KEEP_TEST_TMP:-}" ]]; then echo "kept test directory: $case_dir" >&2; else rm -rf "$case_dir"; fi' EXIT

cat > "$case_dir/fake_pw.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
test "$1" = "-inp"
input="$2"
outdir="$(sed -n "s/^[[:space:]]*outdir[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/ip" "$input")"
prefix="$(sed -n "s/^[[:space:]]*prefix[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/ip" "$input")"
prefix="${prefix:-pwscf}"
save="$outdir/$prefix.save"
mkdir -p "$save"
if grep -qi "startingpot[[:space:]]*=[[:space:]]*['\"]file['\"]" "$input"; then
  test -s "$save/charge-density.dat"
else
  printf 'reference density\n' > "$save/charge-density.dat"
fi
cat <<'OUTPUT'
!    total energy              =     -1.0000000000 Ry
Forces acting on atoms
     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
JOB DONE.
OUTPUT
EOF
chmod +x "$case_dir/fake_pw.sh"

run_case() {
  local name="$1"
  mkdir -p "$case_dir/$name"
  cat > "$case_dir/$name/fdvib.in" <<EOF
scf_input = scf.in
outdir = fdvib
system_type = local
selected_atoms = 1
displacement_angstrom = 0.529177210544
pw_command = bash $case_dir/fake_pw.sh
prefix = system
run_dynmat = false
EOF
  "$fdvib" -inp "$case_dir/$name/fdvib.in" > "$case_dir/$name/run.out"
}

mkdir -p "$case_dir/crystal"
cat > "$case_dir/crystal/scf.in" <<'EOF'
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  prefix = 'u'
  outdir = './out'
/
&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
/
ATOMIC_SPECIES
H 1.0 H.upf
ATOMIC_POSITIONS {crystal}
H 0.25 0.25 0.25
CELL_PARAMETERS (bohr)
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
EOF
run_case crystal
grep -q "H          0.3500000000       0.2500000000       0.2500000000" "$case_dir/crystal/fdvib/calculations/disp_0001_x_p_001/pw.in"
grep -q "H          0.1500000000       0.2500000000       0.2500000000" "$case_dir/crystal/fdvib/calculations/disp_0001_x_m_001/pw.in"

mkdir -p "$case_dir/bohr"
cat > "$case_dir/bohr/scf.in" <<'EOF'
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  prefix = 'u'
  outdir = './out'
/
&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
/
ATOMIC_SPECIES
H 1.0 H.upf
ATOMIC_POSITIONS (bohr)
H 1.0 2.0 3.0
CELL_PARAMETERS {angstrom}
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
EOF
run_case bohr
grep -q "H          2.0000000000       2.0000000000       3.0000000000" "$case_dir/bohr/fdvib/calculations/disp_0001_x_p_001/pw.in"

mkdir -p "$case_dir/alat"
cat > "$case_dir/alat/scf.in" <<'EOF'
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  prefix = 'u'
  outdir = './out'
/
&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
  A = 5.29177210544
/
ATOMIC_SPECIES
H 1.0 H.upf
ATOMIC_POSITIONS {alat}
H 0.1 0.2 0.3
CELL_PARAMETERS alat
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
EOF
run_case alat
grep -q "H          0.2000000000       0.2000000000       0.3000000000" "$case_dir/alat/fdvib/calculations/disp_0001_x_p_001/pw.in"
