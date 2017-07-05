cmake_minimum_required(VERSION 3.9)
project(Polymer)
set(CMAKE_MACOSX_RPATH 1)
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
add_library(atom SHARED src/atom.cpp src/vector.cpp)
add_library(bond SHARED src/bond.cpp src/atom.cpp src/vector.cpp)
add_library(ljpotential SHARED src/vector src/ljpotential.cpp)
# Build executables for testing
add_executable(vector_test src/vector.cpp src/vector_test.cpp)
add_executable(atom_test src/atom.cpp src/atom_test.cpp)
add_executable(bond_test src/bond.cpp src/bond_test.cpp)
add_executable(molecule_test src/molecule.cpp src/molecule_test.cpp)
add_executable(ljpotential_test src/ljpotential.cpp src/ljpotential_test.cpp)
# Link with what we want to test and the GTest and pthread libraries
target_link_libraries(vector_test gtest_main)
target_link_libraries(atom_test gtest_main vector)
target_link_libraries(bond_test gtest_main vector atom)
target_link_libraries(molecule_test gtest_main vector atom bond)
target_link_libraries(ljpotential_test gtest_main vector ljpotential)
