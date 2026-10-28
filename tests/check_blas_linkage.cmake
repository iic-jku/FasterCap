# Fails if a macOS executable links a BLAS other than Accelerate's, whose headers LinAlgebra uses.
# Usage: cmake -DEXECUTABLE=<path> -P check_blas_linkage.cmake

execute_process(COMMAND otool -L ${EXECUTABLE} OUTPUT_VARIABLE libraries RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "otool -L ${EXECUTABLE} failed")
endif()

# without the first line, the executable itself, whose path may contain "blas"
string(FIND "${libraries}" "\n" first_line_end)
math(EXPR first_line_end "${first_line_end} + 1")
string(SUBSTRING "${libraries}" ${first_line_end} -1 libraries)
string(REPLACE "\n" ";" lines "${libraries}")
set(foreign)
foreach(line IN LISTS lines)
    string(TOLOWER "${line}" lower)
    if(lower MATCHES "blas" AND NOT line MATCHES "/System/Library/Frameworks/Accelerate.framework/")
        string(STRIP "${line}" line)
        list(APPEND foreign "${line}")
    endif()
endforeach()
if(foreign)
    message(FATAL_ERROR "${EXECUTABLE} links a BLAS other than Accelerate: ${foreign}")
endif()
