#!/usr/bin/env bash
set -euo pipefail

fdvib="$1"
case_dir="$(mktemp -d)"
trap 'rm -rf "$case_dir"' EXIT

cat > "$case_dir/co.dynG" <<'EOF'
Dynamical matrix file
CO fake Gaussian test
  2    2   0  10.0000000 0 0 0 0 0
Basis vectors
  1 0 0
  0 1 0
  0 0 1
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
append_mode 3 1.0 0.20 0.30
append_mode 4 2.0 0.40 0.50
append_mode 5 3.0 0.60 0.70
append_mode 6 2143.0 -0.80 0.90

"$fdvib" fakeg "$case_dir" > "$case_dir/fakeg.out"
grep -q "Wrote .*co_fake.out" "$case_dir/fakeg.out"
test -s "$case_dir/co_fake.out"
grep -q 'Standard orientation:' "$case_dir/co_fake.out"
grep -q 'Harmonic frequencies' "$case_dir/co_fake.out"
grep -q 'Frequencies --     -1.0000' "$case_dir/co_fake.out"
grep -q '2143.0000' "$case_dir/co_fake.out"
grep -q 'IR Inten    --      0.0000' "$case_dir/co_fake.out"
grep -q ' Atom  AN      X      Y      Z' "$case_dir/co_fake.out"
grep -Eq '^     1   6' "$case_dir/co_fake.out"
grep -Eq '^     2   8' "$case_dir/co_fake.out"
grep -q 'Normal termination of Gaussian' "$case_dir/co_fake.out"

if "$fdvib" fakeg "$case_dir" > /dev/null 2> "$case_dir/error"; then
  echo "fakeg unexpectedly overwrote an existing file" >&2
  exit 1
fi
grep -q 'Refuse to overwrite' "$case_dir/error"
