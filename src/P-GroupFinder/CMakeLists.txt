# src/P-GroupFinder/CMakeLists.txt

include(${CMAKE_CURRENT_SOURCE_DIR}/sources.cmake)

add_executable(P-GroupFinder ${P_GROUPFINDER_SOURCES})

target_compile_definitions(P-GroupFinder PRIVATE ${ENZO_COMPILE_DEFINITIONS})

target_include_directories(P-GroupFinder PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/nrsrc
    ${CMAKE_SOURCE_DIR}/src/enzo
)

if (ENZO_USE_HDF4)
    target_link_libraries(P-GroupFinder PRIVATE HDF4::HDF4)
else()
    target_link_libraries(P-GroupFinder PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(P-GroupFinder PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

if (ENZO_USE_MPI)
    target_link_libraries(P-GroupFinder PRIVATE MPI::MPI_C)
endif()

target_link_libraries(P-GroupFinder PRIVATE m dl)
