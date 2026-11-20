if(FASTFIELDSOLVERS_HEADLESS)
    # Two tiny patches eight units apart cannot be compressed against a patch
    # inside their extent, regardless of the area-equivalent triangle's size.
    # Four physical links must reproduce the direct leaf kernels to roundoff
    # (1e-12 after removing the Coulomb constant), including neutral sources.
    # At separation 64 the same group must remain compressed into two links.
    add_executable(test_conductor_hierarchy test_conductor_hierarchy.cpp)
    target_link_libraries(test_conductor_hierarchy FasterCapSolver)
    add_test(NAME conductor_hierarchy COMMAND test_conductor_hierarchy)
    set_tests_properties(conductor_hierarchy PROPERTIES TIMEOUT 10)
endif()
