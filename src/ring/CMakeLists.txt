# src/ring/CMakeLists.txt

include(${CMAKE_CURRENT_SOURCE_DIR}/sources.cmake)

add_executable(ring ${RING_SOURCES})

target_compile_definitions(ring PRIVATE ${ENZO_COMPILE_DEFINITIONS})

target_include_directories(ring PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/enzo
)

if (ENZO_USE_HDF4)
    target_link_libraries(ring PRIVATE HDF4::HDF4)
else()
    target_link_libraries(ring PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(ring PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

if (ENZO_USE_MPI)
    target_link_libraries(ring PRIVATE MPI::MPI_CXX MPI::MPI_C)
endif()

target_link_libraries(ring PRIVATE m dl)
