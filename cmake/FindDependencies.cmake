# cmake/FindDependencies.cmake
# Dependency management for Enzo

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# HDF5 Dependency
if (ENZO_USE_HDF4)
    find_package(HDF4 REQUIRED)
else()
    find_package(HDF5 REQUIRED COMPONENTS C)
endif()

# MPI Dependency
if (ENZO_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
endif()

# UUID Dependency
if (ENZO_USE_UUID)
    find_package(PkgConfig QUIET)
    if (PKG_CONFIG_FOUND)
        pkg_check_modules(UUID uuid)
    endif()
    if (NOT UUID_FOUND)
        find_path(UUID_INCLUDE_DIR NAMES uuid/uuid.h HINTS /usr/include /usr/local/include)
        find_library(UUID_LIBRARY NAMES uuid HINTS /usr/lib /usr/lib64 /usr/local/lib)
        if (UUID_LIBRARY AND UUID_INCLUDE_DIR)
            set(UUID_FOUND TRUE)
            set(UUID_LIBRARIES ${UUID_LIBRARY})
            set(UUID_INCLUDE_DIRS ${UUID_INCLUDE_DIR})
        endif()
    endif()
endif()

# Grackle Dependency
if (ENZO_USE_GRACKLE)
    find_package(Grackle REQUIRED)
endif()

# HYPRE Dependency
if (ENZO_USE_HYPRE)
    find_package(HYPRE REQUIRED)
endif()

# Python Dependency
if (ENZO_USE_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
endif()
