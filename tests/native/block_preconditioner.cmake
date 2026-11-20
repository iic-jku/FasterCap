if(FASTFIELDSOLVERS_HEADLESS)
    # Integer blocks act exactly on both halves of a complex vector. The
    # material-specific diagonal inverse must compose to identity within 1e-12.
    add_executable(test_block_preconditioner test_block_preconditioner.cpp)
    target_link_libraries(test_block_preconditioner FasterCapSolver)
    add_test(NAME block_preconditioner COMMAND test_block_preconditioner)
    set_tests_properties(block_preconditioner PROPERTIES TIMEOUT 10)
endif()
