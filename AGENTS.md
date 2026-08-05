# Repository Guidelines

## Project Structure & Module Organization

FDVIB is a Go CLI. The entry point is `cmd/fdvib/`; implementation packages
live in `internal/` and mirror the former C++ module layout: `config` and
`settings` (fdvib.in parsing and validation), `qeinput`/`qeoutput`/
`qevibration` (Quantum ESPRESSO parsing), `run` (calculation orchestration
with locking and recovery), `hessian` (force-constant assembly and the .dynG
writer), `modes` (rigid-body classification and mode selection), `export`
(Molden, Shermo, and thermochemistry output), plus `elements`, `process`,
`state`, `results`, `units`, and `init`. The `internal/fakeqe` test doubles
simulate `pw.x` and `dynmat.x` for integration tests. Templates are in
`examples/local/` and `examples/gas/`. GitHub Pages sources are in `docs/`:
keep operational guidance in `index.md`, and equations, algorithms, and
physical limitations in `theory.md` without duplication. Releases use
`.github/workflows/release.yml`.

## Build, Test, and Development Commands

```bash
go build ./...
go vet ./...
go test ./...
go build -ldflags "-X main.version=1.0.4" -o fdvib ./cmd/fdvib
./fdvib --help
```

The first three commands compile, vet, and run the test suite (unit tests plus
end-to-end integration tests driven by `internal/fakeqe`). The fourth builds
the binary with the release version string; run it for CLI checks. Run
`git diff --check` before committing.

## Coding Style & Naming Conventions

Follow surrounding Go style; run `gofmt` and do not reformat unrelated code.
Packages use lowercase names, exported identifiers are `PascalCase`, and
package-private helpers stay lowercase. Keep cross-package data types small
and explicit; return errors with `fmt.Errorf` and keep messages actionable.
The formatting helpers in `internal/config` (FormatGeneral/FormatSci/
FormatFixed) reproduce libstdc++ iostream output byte-for-byte — do not
replace them with `strconv` calls directly. Preserve source SPDX and copyright
headers.

## Testing Guidelines

Code changes require `go vet ./...` and `go test ./...` to pass, plus focused
smoke checks in a temporary directory, including affected error paths and
overwrite protection. Format-sensitive output (dynG, molden, shm, thermo.dat,
forces.dat, state markers) is covered by byte-exact expectations; `testdata/`
holds the C++ iostream formatting reference (`fmt_cpp.txt`). For
documentation-only changes, check links, rendered structure, and
`git diff --check`. Do not commit generated calculation data.

## Commit & Pull Request Guidelines

Use short, lowercase, imperative commit summaries, such as
`align local atom selection templates`; omit version numbers. Keep commits
focused and preserve unrelated worktree changes. Do not push, tag, or rewrite
history unless explicitly requested. Pull requests should explain behavior and
validation. Update README, docs, examples, and changelog where applicable.
Keep the release version, changelog, and release tag aligned.

## Data Safety

Treat the user-selected `scf_input` as read-only. Apply coordinate, path,
`startingpot`, and scratch-directory changes only to generated attempt inputs.
Never overwrite user inputs or completed calculation snapshots silently.
Preserve failure-recovery markers and dataset-digest checks, and keep generated
QE attempts isolated under the configured FDVIB output directory.
