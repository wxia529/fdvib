#!/usr/bin/env bash
set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'rm -rf "$case_dir"' EXIT

cat > "$case_dir/scf.in" <<'EOF'
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  outdir = './out'
/
&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
/
ATOMIC_SPECIES
H 1.00784 H.upf
ATOMIC_POSITIONS angstrom
H 0.0 0.0 0.0
CELL_PARAMETERS angstrom
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
EOF

cat > "$case_dir/fdvib.in" <<'EOF'
scf_input = 'scf.in'
workdir = 'fdvib'
system_type = 'local'
selected_atoms = 1
displacement_angstrom = 0.01
pw_command = 'bash fake_pw.sh'
EOF

cat > "$case_dir/fake_pw.sh" <<'EOF'
#!/usr/bin/env bash
cat <<'OUTPUT'
Forces acting on atoms
     atom    1 type  1   force =     0.10000000   -0.20000000    0.30000000
JOB DONE.
OUTPUT
EOF

"$fdvib" prepare "$case_dir/fdvib.in" > "$case_dir/prepare.out"
grep -q '^Prepared 6 displacement jobs' "$case_dir/prepare.out"
test "$(find "$case_dir/fdvib" -name 'disp_*.in' | wc -l)" -eq 6
test "$(grep -vc '^#' "$case_dir/fdvib/jobs.list")" -eq 6
"$fdvib" prepare "$case_dir/fdvib.in" > /dev/null
"$fdvib" run "$case_dir/fdvib.in" > "$case_dir/run.out"
test "$(find "$case_dir/fdvib" -name forces.dat | wc -l)" -eq 6
grep -q '^1 1.000000000000000e-01 -2.000000000000000e-01 3.000000000000000e-01$' \
  "$case_dir/fdvib/disp_0001_x_p/forces.dat"
"$fdvib" run "$case_dir/fdvib.in" > "$case_dir/restart.out"
grep -q '^Completed 0, preserved 6 existing jobs$' "$case_dir/restart.out"

printf "unknown_parameter = 1\n" >> "$case_dir/fdvib.in"
if "$fdvib" prepare "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "prepare unexpectedly accepted an unknown parameter" >&2
  exit 1
fi
grep -q 'Unknown parameter in fdvib.in' "$case_dir/error"
