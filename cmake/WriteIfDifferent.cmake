# cmake/WriteIfDifferent.cmake
# Write a file only when its content changes, so timestamps (and rebuilds) stay stable.
function(write_if_different PATH CONTENT)
    set(_old "")
    if (EXISTS "${PATH}")
        file(READ "${PATH}" _old)
    endif()
    if (NOT "${_old}" STREQUAL "${CONTENT}")
        file(WRITE "${PATH}" "${CONTENT}")
    endif()
endfunction()
