if(FASTFIELDSOLVERS_HEADLESS)
    # Uniform fields are exact under grouping. Dyadic field components and axial
    # normals allow only roundoff (1e-12), including cancellation of parent normals.
    # A known harmonic quadratic field must integrate exactly after translations.
    # Use exact vertices with rounded float centers, including a 2^20 translation;
    # their polynomial integral remains exact to roundoff (1e-12).
    # For a neutral and nonneutral source pair, halving size/distance must reduce
    # the error by at least 1/0.3: second-order relative convergence gives 1/4,
    # with room for higher-order terms and float centroid coordinates.
    add_executable(test_dielectric_hierarchy test_dielectric_hierarchy.cpp)
    target_link_libraries(test_dielectric_hierarchy FasterCapSolver)
    add_test(NAME dielectric_hierarchy COMMAND test_dielectric_hierarchy)
    set_tests_properties(dielectric_hierarchy PROPERTIES TIMEOUT 10)
endif()
