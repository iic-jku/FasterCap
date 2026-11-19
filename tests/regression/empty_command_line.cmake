# Empty input must print usage and exit 64 before any solver work.
if(FASTFIELDSOLVERS_HEADLESS)
    add_test(NAME empty_command_line
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_empty_command_line.py
            $<TARGET_FILE:FasterCap>)
    set_tests_properties(empty_command_line PROPERTIES TIMEOUT 40)
endif()
