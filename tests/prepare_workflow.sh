#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Wanting Xia

set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'if [[ -n "${KEEP_TEST_TMP:-}" ]]; then echo "kept test directory: $case_dir" >&2; else rm -rf "$case_dir"; fi' EXIT

cat > "$case_dir/scf.in" <<'EOF'
&CONTROL
  calculation = 'scf'
  tprnfor = .true.
  prefix = 'seedtest'
  outdir = './out'
  pseudo_dir = '~/pseudo#hash'
  disk_io = 'high'
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
scf_input = scf.in
outdir = fdvib
system_type = local
selected_atoms = 1
displacement_angstrom = 0.01
pw_command = bash @FAKE_PW@
prefix = system
run_dynmat = false
dynmat_command = bash @FAKE_DYNMAT@
EOF

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
  grep -qx 'reference PAW data' "$save/paw.txt"
  printf 'updated density\n' > "$save/charge-density.dat"
  printf 'updated PAW data\n' > "$save/paw.txt"
else
  printf 'reference density\n' > "$save/charge-density.dat"
  printf 'reference PAW data\n' > "$save/paw.txt"
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
sed -i "s|@FAKE_PW@|$case_dir/fake_pw.sh|" "$case_dir/fdvib.in"

cat > "$case_dir/fake_dynmat.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
test "$1" = "-inp"
filout="$(sed -n "s/^[[:space:]]*filout[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/ip" "$2")"
: > "$filout"
for index in 1 2 3; do
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$index" >> "$filout"
  printf ' ( 1 0 0 0 0 0 )\n' >> "$filout"
done
printf 'JOB DONE.\n'
EOF
sed -i "s|@FAKE_DYNMAT@|$case_dir/fake_dynmat.sh|" "$case_dir/fdvib.in"

if "$fdvib" -in "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted the removed -in option" >&2
  exit 1
fi
grep -q 'fdvib -inp fdvib.in' "$case_dir/error"

cat > "$case_dir/old-style.in" <<'EOF'
&FDVIB
  scf_input = 'scf.in'
/
EOF
if "$fdvib" -inp "$case_dir/old-style.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted namelist fdvib.in syntax" >&2
  exit 1
fi
grep -q 'Namelist syntax is not accepted' "$case_dir/error"

"$fdvib" -inp "$case_dir/fdvib.in" > "$case_dir/run.out"
grep -q '^Running unperturbed reference SCF$' "$case_dir/run.out"
grep -q '^Completed 6, preserved 0 displacement jobs$' "$case_dir/run.out"
grep -q "dynmat.x was not requested" "$case_dir/run.out"
grep -q '^FDVIB calculation completed$' "$case_dir/run.out"
test "$(find "$case_dir/fdvib/calculations" -maxdepth 2 -name pw.in | wc -l)" -eq 6
test "$(grep -Ril "outdir[[:space:]]*=[[:space:]]*'./out'" "$case_dir/fdvib/calculations" --include='*.in' | wc -l)" -eq 7
grep -q "pseudo_dir = '$HOME/pseudo#hash'" "$case_dir/fdvib/calculations/init_scf_001/scf.in"
grep -q "pseudo_dir = '$HOME/pseudo#hash'" "$case_dir/fdvib/calculations/disp_0001_x_p_001/pw.in"
test "$(grep -Ril "startingpot[[:space:]]*=[[:space:]]*'file'" "$case_dir/fdvib/calculations" --include=pw.in | wc -l)" -eq 6
test "$(grep -Ril "disk_io[[:space:]]*=[[:space:]]*'nowf'" "$case_dir/fdvib/calculations" --include=pw.in | wc -l)" -eq 6
grep -q "disk_io = 'high'" "$case_dir/fdvib/calculations/init_scf_001/scf.in"
test "$(find "$case_dir/fdvib/calculations" -path '*/seedtest.save/paw.txt' | wc -l)" -eq 7
test -f "$case_dir/fdvib/state/init_scf.complete"
grep -q '^init_scf_001 charge-density.dat ' "$case_dir/fdvib/state/init_scf.complete"
test -f "$case_dir/fdvib/state/disp_0001_x_m.complete"
if grep -qi "startingpot[[:space:]]*=[[:space:]]*['\"]file['\"]" "$case_dir/fdvib/calculations/init_scf_001/scf.in"; then
  echo "reference SCF unexpectedly uses startingpot='file'" >&2
  exit 1
fi
test -s "$case_dir/fdvib/results/system.dynG"
test -s "$case_dir/fdvib/results/metadata.dat"
grep -q '^program = qe$' "$case_dir/fdvib/results/metadata.dat"
grep -q '^mode_selection = local$' "$case_dir/fdvib/results/metadata.dat"
test ! -e "$case_dir/fdvib/results/fdvib.in.reference"
test ! -e "$case_dir/fdvib/results/electronic_structure.dat"
test ! -e "$case_dir/fdvib/results/system.freq.out"

mkdir "$case_dir/no-disk-io"
cp "$case_dir/scf.in" "$case_dir/no-disk-io/scf.in"
cp "$case_dir/fdvib.in" "$case_dir/no-disk-io/fdvib.in"
sed -i '/disk_io/d' "$case_dir/no-disk-io/scf.in"
sed -i 's|~/pseudo#hash|./pseudo|' "$case_dir/no-disk-io/scf.in"
"$fdvib" -inp "$case_dir/no-disk-io/fdvib.in" > "$case_dir/no-disk-io/run.out"
test "$(grep -Ril "disk_io[[:space:]]*=[[:space:]]*'nowf'" "$case_dir/no-disk-io/fdvib/calculations" --include=pw.in | wc -l)" -eq 6
grep -q "pseudo_dir = '../../../pseudo'" "$case_dir/no-disk-io/fdvib/calculations/disp_0001_x_p_001/pw.in"
if grep -qi 'disk_io' "$case_dir/no-disk-io/fdvib/calculations/init_scf_001/scf.in"; then
  echo "Initial SCF unexpectedly received disk_io" >&2
  exit 1
fi

sed -i 's/run_dynmat = false/run_dynmat = true/' "$case_dir/fdvib.in"
"$fdvib" -inp "$case_dir/fdvib.in" > "$case_dir/restart.out"
grep -q '^Preserved completed reference SCF$' "$case_dir/restart.out"
grep -q '^Completed 0, preserved 6 displacement jobs$' "$case_dir/restart.out"
grep -q '^Preserved completed Hessian analysis$' "$case_dir/restart.out"
grep -q '^Running dynmat.x$' "$case_dir/restart.out"
test -s "$case_dir/fdvib/results/system.freq.out"
test -d "$case_dir/fdvib/calculations/dynmat_001"
grep -q '^dynmat_001 ' "$case_dir/fdvib/state/dynmat.complete"

"$fdvib" -inp "$case_dir/fdvib.in" > "$case_dir/restart-again.out"
grep -q '^Preserved completed dynmat.x result$' "$case_dir/restart-again.out"

rm "$case_dir/fdvib/state/disp_0001_x_m.complete"
rm "$case_dir/fdvib/calculations/disp_0001_x_m_001/forces.dat"
"$fdvib" -inp "$case_dir/fdvib.in" > "$case_dir/recover-displacement.out"
grep -q '^Recovered completed disp_0001_x_m from disp_0001_x_m_001$' "$case_dir/recover-displacement.out"
test ! -e "$case_dir/fdvib/calculations/disp_0001_x_m_002"

rm "$case_dir/fdvib/state/init_scf.complete"
rm "$case_dir/fdvib/metadata.dat"
"$fdvib" -inp "$case_dir/fdvib.in" > "$case_dir/recover-init.out"
grep -q '^Recovered completed initial SCF from init_scf_001$' "$case_dir/recover-init.out"
test ! -e "$case_dir/fdvib/calculations/init_scf_002"

reference_paw="$case_dir/fdvib/calculations/init_scf_001/out/seedtest.save/paw.txt"
printf 'damaged PAW data\n' > "$reference_paw"
if "$fdvib" -inp "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted modified reference PAW data" >&2
  exit 1
fi
grep -q 'Completed reference PAW data is missing or modified' "$case_dir/error"
printf 'reference PAW data\n' > "$reference_paw"

rm "$case_dir/fdvib/state/dynmat.complete"
"$fdvib" -inp "$case_dir/fdvib.in" > "$case_dir/recover.out"
grep -q '^Recovered completed dynmat.x result from dynmat_001$' "$case_dir/recover.out"

"$fdvib" shm "$case_dir/fdvib/results" > "$case_dir/shm.out"
grep -q '^-5.000000000000000E+00$' "$case_dir/fdvib/results/system.shm"
grep -q 'SHM mode selection: local, retained 3, removed 0' "$case_dir/shm.out"

sed -i 's/^format=2$/format=1/' "$case_dir/fdvib/state/dataset.state"
if "$fdvib" -inp "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted the previous dataset layout" >&2
  exit 1
fi
grep -q 'Dataset differs from the existing calculation' "$case_dir/error"
sed -i 's/^format=1$/format=2/' "$case_dir/fdvib/state/dataset.state"

sed -i "/&ELECTRONS/a\\  startingpot = 'file'" "$case_dir/scf.in"
if "$fdvib" -inp "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted startingpot='file' in scf.in" >&2
  exit 1
fi
grep -q "must not set startingpot='file'" "$case_dir/error"

sed -i "/startingpot = 'file'/d" "$case_dir/scf.in"
printf "unknown_parameter = 1\n" >> "$case_dir/fdvib.in"
if "$fdvib" -inp "$case_dir/fdvib.in" > /dev/null 2> "$case_dir/error"; then
  echo "FDVIB unexpectedly accepted an unknown parameter" >&2
  exit 1
fi
grep -q 'Unknown parameter in fdvib.in' "$case_dir/error"
