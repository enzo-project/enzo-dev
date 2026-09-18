# cmake/FindGrackle.cmake
# Find the Grackle chemistry and cooling library

find_path(GRACKLE_INCLUDE_DIR
    NAMES grackle.h
    HINTS ${GRACKLE_ROOT} $ENV{GRACKLE_ROOT} $ENV{HOME}/local $ENV{HOME}/grackle $ENV{HOME}/.local /usr /usr/local
    PATH_SUFFIXES include
)

find_library(GRACKLE_LIBRARY
    NAMES grackle
    HINTS ${GRACKLE_ROOT} $ENV{GRACKLE_ROOT} $ENV{HOME}/local $ENV{HOME}/grackle $ENV{HOME}/.local /usr /usr/local
    PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Grackle
    REQUIRED_VARS GRACKLE_LIBRARY GRACKLE_INCLUDE_DIR
)

if (GRACKLE_FOUND)
    set(GRACKLE_LIBRARIES ${GRACKLE_LIBRARY})
    set(GRACKLE_INCLUDE_DIRS ${GRACKLE_INCLUDE_DIR})
    if (NOT TARGET Grackle::Grackle)
        add_library(Grackle::Grackle UNKNOWN IMPORTED)
        set_target_properties(Grackle::Grackle PROPERTIES
            IMPORTED_LOCATION "${GRACKLE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GRACKLE_INCLUDE_DIR}"
        )
    endif()
endif()
