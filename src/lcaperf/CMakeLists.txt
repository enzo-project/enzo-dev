# src/lcaperf/CMakeLists.txt

add_library(lcaperf_base STATIC jb.C)

set(LCAPERF_SIMPLE_TOOLS
    jb-attribute
    jb-columns
    jb-counter
    jb-delete
    jb-global
    jb-header
    jb-merge
    jb-print
    jb-promote
    jb-region
    jb-rename
)

foreach(tool ${LCAPERF_SIMPLE_TOOLS})
    add_executable(${tool} ${tool}.C)
    target_link_libraries(${tool} PRIVATE lcaperf_base)
endforeach()

add_executable(jb-derive jb-derive.C sort.C)
target_link_libraries(jb-derive PRIVATE lcaperf_base)
