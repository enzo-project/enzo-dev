# src/enzohop/CMakeLists.txt

include(${CMAKE_CURRENT_SOURCE_DIR}/sources.cmake)

add_executable(enzohop ${ENZOHOP_SOURCES})

target_compile_definitions(enzohop PRIVATE ${ENZO_COMPILE_DEFINITIONS})

target_include_directories(enzohop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/enzo
)

if (ENZO_USE_HDF4)
    target_link_libraries(enzohop PRIVATE HDF4::HDF4)
else()
    target_link_libraries(enzohop PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(enzohop PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

if (ENZO_USE_MPI)
    target_link_libraries(enzohop PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
endif()

target_link_libraries(enzohop PRIVATE enzo_lib m dl)
