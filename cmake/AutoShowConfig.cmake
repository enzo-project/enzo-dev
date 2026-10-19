# cmake/AutoShowConfig.cmake
# Function to generate auto_show_*.C files for metadata target

function(generate_auto_show_files TARGET_DIR)
    include(WriteIfDifferent)

    # auto_show_version.C is regenerated at build time (see GenerateVersion.cmake)
    # so the reported git revision tracks commits made after configure.
    add_custom_target(enzo_version
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DOUTPUT_FILE=${TARGET_DIR}/auto_show_version.C
            -DMODULE_DIR=${CMAKE_SOURCE_DIR}/cmake
            -P ${CMAKE_SOURCE_DIR}/cmake/GenerateVersion.cmake
        BYPRODUCTS ${TARGET_DIR}/auto_show_version.C
        COMMENT "Generating auto_show_version.C"
        VERBATIM)
    set_source_files_properties(${TARGET_DIR}/auto_show_version.C PROPERTIES GENERATED TRUE)

    # Generate auto_show_config.C
    write_if_different("${TARGET_DIR}/auto_show_config.C"
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
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _BT)
    write_if_different("${TARGET_DIR}/auto_show_flags.C"
"#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,\"\\n\");
   fprintf (fp,\"CC  = ${CMAKE_C_COMPILER}\\n\");
   fprintf (fp,\"CXX = ${CMAKE_CXX_COMPILER}\\n\");
   fprintf (fp,\"FC  = ${CMAKE_Fortran_COMPILER}\\n\");
   fprintf (fp,\"\\n\");
   fprintf (fp,\"CFLAGS   = ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${_BT}}\\n\");
   fprintf (fp,\"CXXFLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_BT}}\\n\");
   fprintf (fp,\"FFLAGS   = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_BT}}\\n\");
   fprintf (fp,\"\\n\");
}
"
    )
endfunction()
