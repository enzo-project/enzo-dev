# Root CMakeLists.txt for Enzo
cmake_minimum_required(VERSION 3.20)

project(Enzo LANGUAGES C CXX Fortran)

# Path for custom CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set build output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# Configure compiler flags for Fortran
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    include(CheckFortranCompilerFlag)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-second-underscore -ffixed-line-length-132 -std=legacy")
    check_fortran_compiler_flag("-fallow-argument-mismatch" HAS_ALLOW_ARGUMENT_MISMATCH)
    if (HAS_ALLOW_ARGUMENT_MISMATCH)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-literal-suffix -Wno-write-strings")
endif()

# Include configuration options and preprocessor definitions
include(EnzoOptions)

# Include dependency search
include(FindDependencies)

# Include metadata generator
include(AutoShowConfig)

# Add source directory
add_subdirectory(src)
