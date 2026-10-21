# cmake/GenerateVersion.cmake
# Run at BUILD time (cmake -P) to regenerate auto_show_version.C from git.
# Expects: SOURCE_DIR, OUTPUT_FILE, MODULE_DIR
include("${MODULE_DIR}/WriteIfDifferent.cmake")

set(GIT_BRANCH "")
set(GIT_REVISION "")
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${SOURCE_DIR}/.git")
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${SOURCE_DIR}" OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY "${SOURCE_DIR}" OUTPUT_VARIABLE GIT_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
endif()
if (NOT GIT_BRANCH)
    set(GIT_BRANCH "unknown")
endif()
if (NOT GIT_REVISION)
    set(GIT_REVISION "unknown")
endif()

write_if_different("${OUTPUT_FILE}"
"#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,\"\\n\");
   fprintf (fp,\"Git Branch   ${GIT_BRANCH}\\n\");
   fprintf (fp,\"Git Revision ${GIT_REVISION}\\n\");
   fprintf (fp,\"\\n\");
}
")
