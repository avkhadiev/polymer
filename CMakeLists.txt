cmake_minimum_required(VERSION 3.9)
project(Polymer)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build)
include_directories("${gtest_SOURCE_DIR}/include")

# Add header files from include/ directory
include_directories(include)
add_library(vector_and_parsing STATIC src/vector.cpp src/parsing.cpp)
add_library(atom_and_bond STATIC src/vector.cpp src/parsing.cpp src/simple_atom.cpp src/simple_bond.cpp)
add_library(polymer STATIC src/vector.cpp src/parsing.cpp src/simple_atom.cpp src/simple_bond.cpp src/simple_polymer.cpp)
add_library(state STATIC src/vector.cpp src/parsing.cpp src/simple_atom.cpp src/simple_bond.cpp src/simple_polymer.cpp src/simple_state.cpp)
add_library(MD_simulation STATIC src/parsing.cpp src/vector.cpp  src/simple_atom.cpp src/simple_polymer.cpp src/simple_bond.cpp src/simple_state.cpp src/config_handler.cpp src/simple_simulation.cpp  src/observables.cpp src/observable_container src/ljpotential.cpp  src/force_updater.cpp src/integrator.cpp src/verlet_integrator.cpp src/rattle_integrator.cpp)
# Build executables for testing
add_executable(vector_test src/vector_test.cpp src/vector.cpp)
add_executable(parsing_test src/parsing_test.cpp src/parsing.cpp)
add_executable(simple_atom_test src/simple_atom_test.cpp)
add_executable(simple_bond_test src/simple_bond_test.cpp)
add_executable(simple_polymer_test src/simple_polymer_test.cpp)
add_executable(simple_state_test src/simple_state_test.cpp)
add_executable(ljpotential_test src/ljpotential_test.cpp src/ljpotential.cpp)
add_executable(observables_test src/observables_test.cpp src/observables.cpp)
add_executable(observable_container_test src/observable_container_test.cpp src/observable_container.cpp src/observables.cpp)
add_executable(triatomic src/triatomic.cpp)
# Link with what we want to test and the GTest and pthread libraries
target_link_libraries(vector_test gtest_main)
target_link_libraries(parsing_test gtest_main)
target_link_libraries(simple_atom_test gtest_main atom_and_bond)
target_link_libraries(simple_bond_test gtest_main atom_and_bond)
target_link_libraries(simple_polymer_test gtest_main polymer)
target_link_libraries(simple_state_test gtest_main state)
target_link_libraries(ljpotential_test gtest_main vector_and_parsing)
target_link_libraries(observables_test gtest_main vector_and_parsing)
target_link_libraries(observable_container_test gtest_main vector_and_parsing)
target_link_libraries(triatomic gtest_main MD_simulation)