# Repository Guidelines

## Project Structure & Module Organization

FDVIB is a C++17 CLI. Code lives in `src/`; `fdvib.hpp` contains shared
interfaces and focused modules handle QE parsing, execution, force constants,
mode export, and thermochemistry. Templates are in `examples/local/` and
`examples/gas/`. GitHub Pages sources are in `docs/`: keep operational guidance
in `index.md`, and equations, algorithms, and physical limitations in
`theory.md` without duplication. Releases use `.github/workflows/release.yml`.

## Build, Test, and Development Commands

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/fdvib --help
cmake --install build --prefix /tmp/fdvib-install
```

The first two commands configure and compile with warnings enabled. Run the
binary for CLI checks; the install command verifies packaging. Run
`git diff --check` before committing.

## Coding Style & Naming Conventions

Follow surrounding C++ style; do not reformat unrelated code. Declare
cross-module interfaces in `src/fdvib.hpp` and keep private helpers in anonymous
namespaces. Files and functions use `snake_case`, types `PascalCase`, and
constants `UPPER_SNAKE_CASE`. Prefer `std::filesystem`; use POSIX APIs only for
process control, locking, or similar system features. Report actionable errors
with `std::runtime_error`. Preserve source SPDX and copyright headers. No
formatter is configured; compile with `-Wall -Wextra -Wpedantic`.

## Testing Guidelines

There is no automated test target or coverage threshold. Code changes require
a fresh Release build and focused smoke checks in a temporary directory,
including affected error paths and overwrite protection. For documentation-only
changes, check links, rendered structure, and `git diff --check`. Do not commit
generated calculation data.

## Commit & Pull Request Guidelines

Use short, lowercase, imperative commit summaries, such as
`align local atom selection templates`; omit version numbers. Keep commits
focused and preserve unrelated worktree changes. Do not push, tag, or rewrite
history unless explicitly requested. Pull requests should explain behavior and
validation. Update README, docs, examples, and changelog where applicable. Keep
the CMake version, changelog, and release tag aligned.

## Data Safety

Treat the user-selected `scf_input` as read-only. Apply coordinate, path,
`startingpot`, and scratch-directory changes only to generated attempt inputs.
Never overwrite user inputs or completed calculation snapshots silently.
Preserve failure-recovery markers and dataset-digest checks, and keep generated
QE attempts isolated under the configured FDVIB output directory.
