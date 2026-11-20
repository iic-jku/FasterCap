if(FASTFIELDSOLVERS_HEADLESS)
    # The worker path must preserve finite values exactly and reject failed
    # coefficients without publishing a partial/nonfinite value or a counter write.
    add_executable(test_parallel_coefficients test_parallel_coefficients.cpp)
    target_link_libraries(test_parallel_coefficients FasterCapSolver)
    add_test(NAME parallel_coefficients COMMAND test_parallel_coefficients)
    set_tests_properties(parallel_coefficients PROPERTIES TIMEOUT 10)
endif()
