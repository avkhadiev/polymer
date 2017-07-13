cmake_minimum_required(VERSION 3.9)
project(Polymer)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
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
add_library(vector SHARED src/vector.cpp)
add_library(parsing SHARED src/parsing.cpp)
add_library(state SHARED src/parsing.cpp src/vector.cpp src/atom.cpp src/molecule.cpp src/bond.cpp src/state.cpp)
add_library(observable_containers SHARED src/parsing.cpp src/vector.cpp src/observables.cpp src/observable_container)
# Build executables for testing
add_executable(vector_test src/vector_test.cpp)
add_executable(parsing_test src/parsing_test.cpp)
add_executable(atom_test src/atom_test.cpp)
add_executable(bond_test src/bond_test.cpp)
add_executable(molecule_test src/molecule_test.cpp)
add_executable(state_test src/state_test.cpp)
add_executable(ljpotential_test src/ljpotential.cpp src/ljpotential_test.cpp)
add_executable(observables_test src/observables_test.cpp)
add_executable(observable_container_test src/observable_container_test.cpp)
add_executable(initialize_triatomic src/initialize_triatomic.cpp)
# Link with what we want to test and the GTest and pthread libraries
target_link_libraries(vector_test gtest_main vector)
target_link_libraries(parsing_test gtest_main parsing)
target_link_libraries(atom_test gtest_main state)
target_link_libraries(bond_test gtest_main state)
target_link_libraries(molecule_test gtest_main state)
target_link_libraries(state_test gtest_main state)
target_link_libraries(ljpotential_test gtest_main vector parsing)
target_link_libraries(observables_test gtest_main observable_containers)
target_link_libraries(observable_container_test gtest_main observable_containers)
target_link_libraries(initialize_triatomic state)
