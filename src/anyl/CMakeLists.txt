# src/anyl/CMakeLists.txt

include(${CMAKE_CURRENT_SOURCE_DIR}/sources.cmake)

add_executable(anyl ${ANYL_SOURCES})

target_compile_definitions(anyl PRIVATE ${ENZO_COMPILE_DEFINITIONS})

target_include_directories(anyl PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/enzo
)

target_link_libraries(anyl PRIVATE enzo_lib)
