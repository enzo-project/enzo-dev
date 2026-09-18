# src/enzo/CMakeLists.txt

include(${CMAKE_CURRENT_SOURCE_DIR}/sources.cmake)

# Generate metadata files (auto_show_*.C)
generate_auto_show_files("${CMAKE_CURRENT_BINARY_DIR}")

# Include generated auto_show_*.C
set(ENZO_AUTOGEN_SOURCES
    "${CMAKE_CURRENT_BINARY_DIR}/auto_show_config.C"
    "${CMAKE_CURRENT_BINARY_DIR}/auto_show_flags.C"
    "${CMAKE_CURRENT_BINARY_DIR}/auto_show_version.C"
)

set(ENZO_ALL_LIB_SOURCES ${ENZO_CORE_SOURCES})
if (ENZO_USE_CUDA)
    list(APPEND ENZO_ALL_LIB_SOURCES ${ENZO_CUDA_SOURCES})
endif()

# Build enzo core library for reuse by tools (enzohop, anyl, etc.)
add_library(enzo_lib STATIC ${ENZO_ALL_LIB_SOURCES})

target_compile_definitions(enzo_lib PUBLIC ${ENZO_COMPILE_DEFINITIONS})

target_include_directories(enzo_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Compiler flags for Fortran
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(BlockSolve.F FOF_InitialGuess.F90 FSProb_InitialGuess.F90 FSProb_RadiationSource.F90 FSProb_SetupSystem.F90
        PROPERTIES COMPILE_FLAGS "-fno-second-underscore -std=legacy")
endif()

# Dependencies & Linking for enzo_lib
if (ENZO_USE_HDF4)
    target_link_libraries(enzo_lib PUBLIC HDF4::HDF4)
else()
    target_link_libraries(enzo_lib PUBLIC ${HDF5_LIBRARIES})
    target_include_directories(enzo_lib PUBLIC ${HDF5_INCLUDE_DIRS})
endif()

if (ENZO_USE_MPI)
    target_link_libraries(enzo_lib PUBLIC MPI::MPI_CXX MPI::MPI_Fortran)
endif()

if (ENZO_USE_GRACKLE AND TARGET Grackle::Grackle)
    target_link_libraries(enzo_lib PUBLIC Grackle::Grackle)
elseif (ENZO_USE_GRACKLE)
    target_link_libraries(enzo_lib PUBLIC ${GRACKLE_LIBRARIES})
    target_include_directories(enzo_lib PUBLIC ${GRACKLE_INCLUDE_DIRS})
endif()

if (ENZO_USE_HYPRE AND TARGET HYPRE::HYPRE)
    target_link_libraries(enzo_lib PUBLIC HYPRE::HYPRE)
elseif (ENZO_USE_HYPRE)
    target_link_libraries(enzo_lib PUBLIC ${HYPRE_LIBRARIES})
    target_include_directories(enzo_lib PUBLIC ${HYPRE_INCLUDE_DIRS})
endif()

if (ENZO_USE_UUID)
    target_link_libraries(enzo_lib PUBLIC ${UUID_LIBRARIES})
    if (UUID_INCLUDE_DIRS)
        target_include_directories(enzo_lib PUBLIC ${UUID_INCLUDE_DIRS})
    endif()
endif()

target_link_libraries(enzo_lib PUBLIC m dl)

# Build main enzo executable
add_executable(enzo enzo.C ${ENZO_AUTOGEN_SOURCES})
target_link_libraries(enzo PRIVATE enzo_lib)
