# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Wanting Xia

cmake_minimum_required(VERSION 3.12)
project(fdvib VERSION 0.4.3 LANGUAGES CXX)

add_executable(fdvib
  src/main.cpp
  src/common.cpp
  src/qe_input.cpp
  src/qe_output.cpp
  src/qe_vibration.cpp
  src/elements.cpp
  src/process.cpp
  src/prepare_run.cpp
  src/run_state.cpp
  src/analysis.cpp
  src/molden.cpp
  src/thermo.cpp
  src/shm.cpp
  src/mode_selection.cpp
  src/results.cpp
)
target_compile_features(fdvib PRIVATE cxx_std_17)
target_compile_options(fdvib PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(fdvib PRIVATE FDVIB_VERSION="${PROJECT_VERSION}")
install(TARGETS fdvib RUNTIME DESTINATION bin)
install(FILES README.md CHANGELOG.md LICENSE DESTINATION share/doc/fdvib)
install(DIRECTORY docs/ DESTINATION share/doc/fdvib/reference)
install(DIRECTORY examples/ DESTINATION share/fdvib/examples)
