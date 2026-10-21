# src/inits/CMakeLists.txt

include(${CMAKE_CURRENT_SOURCE_DIR}/sources.cmake)

add_executable(inits ${INITS_SOURCES})

target_compile_definitions(inits PRIVATE ${ENZO_COMPILE_DEFINITIONS})

target_include_directories(inits PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/enzo
)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(select_fft.F90 fft90.F90
        PROPERTIES COMPILE_FLAGS "-fno-second-underscore -std=legacy")
endif()

if (ENZO_USE_HDF4)
    target_link_libraries(inits PRIVATE HDF4::HDF4)
else()
    target_link_libraries(inits PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(inits PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

if (ENZO_USE_MPI)
    target_link_libraries(inits PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
endif()

target_link_libraries(inits PRIVATE m dl)
