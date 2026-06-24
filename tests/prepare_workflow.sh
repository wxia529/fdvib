#!/usr/bin/env bash
set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'rm -rf "$case_dir"' EXIT

cat > "$case_dir/scf.in" <<'EOF'
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  prefix = 'seedtest'
  outdir = './out'
/
&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
/
&ELECTRONS
  conv_thr = 1.0d-10
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
&FDVIB
  scf_input = 'scf.in'
  outdir = 'fdvib'
  system_type = 'local'
  selected_atoms = 1
  displacement_angstrom = 0.01
  pw_command = 'bash fake_pw.sh'
  prefix = 'system'
  run_dynmat = .false.
  dynmat_command = 'bash @FAKE_DYNMAT@'
/
EOF

cat > "$case_dir/fake_pw.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
test "$1" = "-in"
input="$2"
outdir="$(sed -n "s/^[[:space:]]*outdir[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/ip" "$input")"
prefix="$(sed -n "s/^[[:space:]]*prefix[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/ip" "$input")"
prefix="${prefix:-pwscf}"
save="$outdir/$prefix.save"
mkdir -p "$save"
if grep -qi "startingpot[[:space:]]*=[[:space:]]*['\"]file['\"]" "$input"; then
  test -s "$save/charge-density.dat"
  printf 'updated density\n' > "$save/charge-density.dat"
else
  printf 'reference density\n' > "$save/charge-density.dat"
fi
cat <<'OUTPUT'
!    total energy              =     -10.0000000000 Ry
Forces acting on atoms
     atom    1 type  1   force =     0.10000000   -0.20000000    0.30000000
The non-local contrib. to forces
     atom    1 type  1   force =     9.00000000    8.00000000    7.00000000
JOB DONE.
OUTPUT
EOF

cat > "$case_dir/fake_dynmat.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
test "$1" = "-in"
filout="$(sed -n "s/^[[:space:]]*filout[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/ip" "$2")"
: > "$filout"
for index in 1 2 3; do
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$index" >> "$filout"
  printf ' ( 1 0 0 0 0 0 )\n' >> "$filout"
done
printf 'JOB DONE.\n'
EOF
sed -i "s|@FAKE_DYNMAT@|$case_dir/fake_dynmat.sh|" "$case_dir/fdvib.in"

"$fdvib" -in "$case_dir/fdvib.in" > "$case_dir/run.out"
grep -q '^Running unperturbed reference SCF$' "$case_dir/run.out"
grep -q '^Completed 6, preserved 0 displacement jobs$' "$case_dir/run.out"
grep -q "dynmat.x was not requested" "$case_dir/run.out"
test "$(find "$case_dir/fdvib/jobs" -name pw.in | wc -l)" -eq 6
test "$(grep -Ril "startingpot[[:space:]]*=[[:space:]]*'file'" "$case_dir/fdvib/jobs" --include=pw.in | wc -l)" -eq 6
if grep -qi "startingpot[[:space:]]*=[[:space:]]*['\"]file['\"]" "$case_dir/fdvib/reference/attempt-001/scf.in"; then
  echo "reference SCF unexpectedly uses startingpot='file'" >&2
  exit 1
fi
test -s "$case_dir/fdvib/results/system.dynG"
test ! -e "$case_dir/fdvib/results/system.freq.out"

sed -i 's/run_dynmat = .false./run_dynmat = .true./' "$case_dir/fdvib.in"
"$fdvib" -in "$case_dir/fdvib.in" > "$case_dir/restart.out"
grep -q '^Preserved completed reference SCF$' "$case_dir/restart.out"
grep -q '^Completed 0, preserved 6 displacement jobs$' "$case_dir/restart.out"
grep -q '^Preserved completed Hessian analysis$' "$case_dir/restart.out"
grep -q '^Running dynmat.x$' "$case_dir/restart.out"
test -s "$case_dir/fdvib/results/system.freq.out"

"$fdvib" -in "$case_dir/fdvib.in" > "$case_dir/restart-again.out"
grep -q '^Preserved completed dynmat.x result$' "$case_dir/restart-again.out"

rm "$case_dir/fdvib/state/dynmat.complete"
"$fdvib" -in "$case_dir/fdvib.in" > "$case_dir/recover.out"
grep -q '^Recovered completed dynmat.x result$' "$case_dir/recover.out"

"$fdvib" shm "$case_dir/fdvib/results" > "$case_dir/shm.out"
grep -q '^-5.000000000000000E+00$' "$case_dir/fdvib/results/system.shm"
grep -q 'SHM mode selection: local, retained 3, removed 0' "$case_dir/shm.out"

sed -i "/&ELECTRONS/a\\  startingpot = 'file'" "$case_dir/scf.in"
if "$fdvib" -in "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted startingpot='file' in scf.in" >&2
  exit 1
fi
grep -q "must not set startingpot='file'" "$case_dir/error"

sed -i "/startingpot = 'file'/d" "$case_dir/scf.in"
printf "unknown_parameter = 1\n" >> "$case_dir/fdvib.in"
if "$fdvib" -in "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted an unknown parameter" >&2
  exit 1
fi
grep -q 'Unknown parameter in fdvib.in' "$case_dir/error"
