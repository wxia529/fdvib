#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Wanting Xia

set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'rm -rf "$case_dir"' EXIT

cat > "$case_dir/co.dynG" <<'EOF'
Dynamical matrix file
CO Molden export test
  2    2   0  10.0000000 0 0 0 0 0
Basis vectors
  1.0 0.0 0.0
  0.5 1.0 0.0
  0.0 0.0 2.0
  1 'C1     ' 10937.3309176688
  2 'O1     ' 14582.1922119821
    1    1 0.0000000000 0.0000000000 0.0000000000
    2    2 0.2000000000 0.0000000000 0.0000000000
EOF

append_mode() {
  local index="$1" frequency="$2" cx="$3" ox="$4"
  printf ' freq ( %s) = 0.0 [THz] = %s [cm-1]\n' "$index" "$frequency" >> "$case_dir/co.freq.out"
  printf ' ( %s 0 0 0 0 0 )\n ( %s 0 0 0 0 0 )\n' "$cx" "$ox" >> "$case_dir/co.freq.out"
}

: > "$case_dir/co.freq.out"
append_mode 1 -1.0 0.10 -0.10
append_mode 2 0.0 0.00 0.00
append_mode 3 0.00000001 0.20 0.30
append_mode 4 2.0 0.40 0.50
append_mode 5 3.0 0.60 0.70
append_mode 6 2143.0 -0.80 0.90

(cd "$case_dir" && "$fdvib" modes . > modes.out)
grep -qx 'Wrote 5 Molden modes to co.mol' "$case_dir/modes.out"
test -s "$case_dir/co.mol"
test ! -e "$case_dir/co.mold"
test ! -e "$case_dir/co_fake.out"

test "$(grep -nE '^\[(Molden Format|Cell|Atoms|FREQ|FR-COORD|FR-NORM-COORD)' "$case_dir/co.mol" | cut -d: -f2- | paste -sd ' ' -)" = \
  '[Molden Format] [Cell] [Atoms] AU [FREQ] [FR-COORD] [FR-NORM-COORD]'
if grep -q '^\[INT\]' "$case_dir/co.mol"; then
  echo "Molden export unexpectedly contains IR intensities" >&2
  exit 1
fi
grep -Eq '^[[:space:]]*5\.29177211[[:space:]]+0\.00000000[[:space:]]+0\.00000000$' "$case_dir/co.mol"
grep -Eq '^[[:space:]]*2\.64588605[[:space:]]+5\.29177211[[:space:]]+0\.00000000$' "$case_dir/co.mol"
grep -Eq '^[[:space:]]*0\.00000000[[:space:]]+0\.00000000[[:space:]]+10\.58354421$' "$case_dir/co.mol"
grep -Eq '^[[:space:]]*C[[:space:]]+1[[:space:]]+6[[:space:]]+0\.00000000' "$case_dir/co.mol"
grep -Eq '^[[:space:]]*O[[:space:]]+2[[:space:]]+8[[:space:]]+2\.00000000' "$case_dir/co.mol"

awk '/^\[FREQ\]/{inside=1;next}/^\[/{inside=0}inside{print}' "$case_dir/co.mol" > "$case_dir/frequencies"
test "$(wc -l < "$case_dir/frequencies")" -eq 5
grep -qx -- '-1.00000000' "$case_dir/frequencies"
grep -qx '0.00000001' "$case_dir/frequencies"
grep -qx '2143.00000000' "$case_dir/frequencies"
if grep -qx '0.00000000' "$case_dir/frequencies"; then
  echo "Molden export unexpectedly kept an exact zero frequency" >&2
  exit 1
fi
test "$(grep -c '^ vibration ' "$case_dir/co.mol")" -eq 5

"$fdvib" modes "$case_dir" > "$case_dir/modes-rerun.out"
grep -q 'Wrote 5 Molden modes to .*co.mol' "$case_dir/modes-rerun.out"

gas_dir="$case_dir/gas-co"
mkdir -p "$gas_dir"
cp "$case_dir/co.dynG" "$gas_dir/co.dynG"
cp "$case_dir/co.freq.out" "$gas_dir/co.freq.out"
cat > "$gas_dir/metadata.dat" <<'EOF'
program = qe
mode_selection = gas
EOF
"$fdvib" modes "$gas_dir" > "$gas_dir/modes.out"
test -s "$gas_dir/co.mol"
if grep -q '^\[Cell\]' "$gas_dir/co.mol"; then
  echo "Gas Molden export unexpectedly contains a cell" >&2
  exit 1
fi
awk '/^\[FREQ\]/{inside=1;next}/^\[/{inside=0}inside{print}' "$gas_dir/co.mol" > "$gas_dir/frequencies"
test "$(wc -l < "$gas_dir/frequencies")" -eq 1
grep -qx '2143.00000000' "$gas_dir/frequencies"

if "$fdvib" fakeg "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "Removed fakeg command unexpectedly succeeded" >&2
  exit 1
fi
if grep -q 'fdvib fakeg' "$case_dir/error"; then
  echo "Usage still advertises removed fakeg command" >&2
  exit 1
fi
