# src/CMakeLists.txt

add_subdirectory(enzo)
add_subdirectory(inits)
add_subdirectory(enzohop)
# ring and P-GroupFinder include mpi.h unconditionally
if (ENZO_USE_MPI)
    add_subdirectory(ring)
endif()
add_subdirectory(anyl)
if (ENZO_USE_MPI)
    add_subdirectory(P-GroupFinder)
endif()
add_subdirectory(lcaperf)
