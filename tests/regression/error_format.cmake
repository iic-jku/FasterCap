# #43: command-line error text is data, including percent conversions.
add_test(NAME error_format
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_error_format.py
        $<TARGET_FILE:FasterCap>)
set_tests_properties(error_format PROPERTIES TIMEOUT 30)
