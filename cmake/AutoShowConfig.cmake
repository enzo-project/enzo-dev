# cmake/AutoShowConfig.cmake
# Function to generate auto_show_*.C files for metadata target

function(generate_auto_show_files TARGET_DIR)
    find_package(Git QUIET)
    if (GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            OUTPUT_VARIABLE GIT_BRANCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            OUTPUT_VARIABLE GIT_REVISION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()

    if (NOT GIT_BRANCH)
        set(GIT_BRANCH "main")
    endif()
    if (NOT GIT_REVISION)
        set(GIT_REVISION "unknown")
    endif()

    # Generate auto_show_version.C
    file(WRITE "${TARGET_DIR}/auto_show_version.C"
"#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,\"\\n\");
   fprintf (fp,\"Git Branch   ${GIT_BRANCH}\\n\");
   fprintf (fp,\"Git Revision ${GIT_REVISION}\\n\");
   fprintf (fp,\"\\n\");
}
"
    )

    # Generate auto_show_config.C
    file(WRITE "${TARGET_DIR}/auto_show_config.C"
"#include <stdio.h>
void auto_show_config(FILE *fp) {
   fprintf (fp,\"\\n\");
   fprintf (fp,\"   MACHINE: CMake Generated Build\\n\");
   fprintf (fp,\"   MACHINE-NAME: cmake\\n\");
   fprintf (fp,\"\\n\");
   fprintf (fp,\"   PARAMETER_MAX_SUBGRIDS  [max-subgrids-###]                : ${ENZO_MAX_SUBGRIDS}\\n\");
   fprintf (fp,\"   PARAMETER_MAX_BARYONS  [max-baryons-###]                  : ${ENZO_MAX_BARYONS}\\n\");
   fprintf (fp,\"   PARAMETER_MAX_TASKS_PER_NODE  [max-tasks-per-node-###]    : ${ENZO_MAX_TASKS_PER_NODE}\\n\");
   fprintf (fp,\"   PARAMETER_MEMORY_POOL_SIZE  [memory-pool-###]             : ${ENZO_MEMORY_POOL_SIZE}\\n\");
   fprintf (fp,\"\\n\");
   fprintf (fp,\"   CONFIG_PRECISION  [precision-{32,64}]                     : ${ENZO_PRECISION}\\n\");
   fprintf (fp,\"   CONFIG_PARTICLES  [particles-{32,64,128}]                 : ${ENZO_PARTICLES}\\n\");
   fprintf (fp,\"   CONFIG_INTEGERS  [integers-{32,64}]                       : ${ENZO_INTEGERS}\\n\");
   fprintf (fp,\"   CONFIG_PARTICLE_IDS  [particle-id-{32,64}]                : ${ENZO_PARTICLE_IDS}\\n\");
   fprintf (fp,\"   CONFIG_INITS  [inits-{32,64}]                             : ${ENZO_INITS_PRECISION}\\n\");
   fprintf (fp,\"   CONFIG_IO  [io-{32,64}]                                   : ${ENZO_IO_PRECISION}\\n\");
   fprintf (fp,\"\\n\");
}
"
    )

    # Generate auto_show_flags.C
    file(WRITE "${TARGET_DIR}/auto_show_flags.C"
"#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,\"\\n\");
   fprintf (fp,\"CC  = ${CMAKE_C_COMPILER}\\n\");
   fprintf (fp,\"CXX = ${CMAKE_CXX_COMPILER}\\n\");
   fprintf (fp,\"FC  = ${CMAKE_Fortran_COMPILER}\\n\");
   fprintf (fp,\"\\n\");
   fprintf (fp,\"CFLAGS   = ${CMAKE_C_FLAGS}\\n\");
   fprintf (fp,\"CXXFLAGS = ${CMAKE_CXX_FLAGS}\\n\");
   fprintf (fp,\"FFLAGS   = ${CMAKE_Fortran_FLAGS}\\n\");
   fprintf (fp,\"\\n\");
}
"
    )
endfunction()
