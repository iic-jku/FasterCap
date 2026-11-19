if(FASTFIELDSOLVERS_HEADLESS)
    # Identical centroids and tied subgroups must produce bounded trees with every leaf once.
    # The 1e-12 area/length allowance checks conservation while adding double-precision values.
    add_executable(test_hierarchy_partition test_hierarchy_partition.cpp)
    target_link_libraries(test_hierarchy_partition FasterCapSolver)
    add_test(NAME hierarchy_partition COMMAND test_hierarchy_partition)
    set_tests_properties(hierarchy_partition PROPERTIES TIMEOUT 10)
endif()
