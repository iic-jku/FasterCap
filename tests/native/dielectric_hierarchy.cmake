if(FASTFIELDSOLVERS_HEADLESS)
    # Uniform fields are exact under grouping. Dyadic field components and axial
    # normals allow only roundoff (1e-12), including cancellation of parent normals.
    add_executable(test_dielectric_hierarchy test_dielectric_hierarchy.cpp)
    target_link_libraries(test_dielectric_hierarchy FasterCapSolver)
    add_test(NAME dielectric_hierarchy COMMAND test_dielectric_hierarchy)
    set_tests_properties(dielectric_hierarchy PROPERTIES TIMEOUT 10)
endif()
