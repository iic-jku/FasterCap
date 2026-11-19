# #44: implicit Jacobi is replaced by an explicit selection; block combinations fail.
add_test(NAME preconditioner_options
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_preconditioner_options.py
        $<TARGET_FILE:FasterCap>)
set_tests_properties(preconditioner_options PROPERTIES TIMEOUT 60)
