cmake_minimum_required(VERSION 3.12)
project(fdvib VERSION 0.4.1 LANGUAGES CXX)

add_executable(fdvib
  src/main.cpp
  src/common.cpp
  src/qe_io.cpp
  src/prepare_run.cpp
  src/analysis.cpp
  src/thermo.cpp
  src/shm.cpp
  src/fakeg.cpp
  src/results.cpp
)
target_compile_features(fdvib PRIVATE cxx_std_17)
target_compile_options(fdvib PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(fdvib PRIVATE FDVIB_VERSION="${PROJECT_VERSION}")
install(TARGETS fdvib RUNTIME DESTINATION bin)
install(FILES README.md CHANGELOG.md LICENSE DESTINATION share/doc/fdvib)
install(DIRECTORY docs/ DESTINATION share/doc/fdvib/reference)
install(DIRECTORY examples/ DESTINATION share/fdvib/examples)

include(CTest)
if(BUILD_TESTING)
  add_test(NAME cli_version COMMAND fdvib --version)
  set_tests_properties(cli_version PROPERTIES PASS_REGULAR_EXPRESSION "fdvib ${PROJECT_VERSION}")
  add_test(NAME gas_rrho_modes
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/gas_rrho_modes.sh $<TARGET_FILE:fdvib>)
  add_test(NAME prepare_workflow
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/prepare_workflow.sh $<TARGET_FILE:fdvib>)
  add_test(NAME qe_units
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/qe_units.sh $<TARGET_FILE:fdvib>)
  add_test(NAME shm_export
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/shm_export.sh $<TARGET_FILE:fdvib>)
  add_test(NAME fakeg_export
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/fakeg_export.sh $<TARGET_FILE:fdvib>)
endif()
