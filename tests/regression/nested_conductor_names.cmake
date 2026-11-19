# #46: small nested inputs must reject oversized decorated names before copying.
add_test(NAME nested_conductor_names
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_nested_names.py
        $<TARGET_FILE:FasterCap>)
set_tests_properties(nested_conductor_names PROPERTIES TIMEOUT 60)
