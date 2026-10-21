# cmake/FindHYPRE.cmake
# Find HYPRE library

find_path(HYPRE_INCLUDE_DIR
    NAMES HYPRE.h
    HINTS ${HYPRE_ROOT} $ENV{HYPRE_ROOT} $ENV{HOME}/local $ENV{HOME}/.local /usr /usr/local
    PATH_SUFFIXES include
)

find_library(HYPRE_LIBRARY
    NAMES HYPRE
    HINTS ${HYPRE_ROOT} $ENV{HYPRE_ROOT} $ENV{HOME}/local $ENV{HOME}/.local /usr /usr/local
    PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE
    REQUIRED_VARS HYPRE_LIBRARY HYPRE_INCLUDE_DIR
)

if (HYPRE_FOUND)
    set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
    set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
    if (NOT TARGET HYPRE::HYPRE)
        add_library(HYPRE::HYPRE UNKNOWN IMPORTED)
        set_target_properties(HYPRE::HYPRE PROPERTIES
            IMPORTED_LOCATION "${HYPRE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIR}"
        )
    endif()
endif()
