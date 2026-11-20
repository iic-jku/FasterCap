if(FASTFIELDSOLVERS_HEADLESS)
    # Known linear systems and independent residuals: tolerance 1e-12 leaves room for
    # roundoff in these small well-conditioned systems, including flexible preconditioning.
    # Automatic residual correction must recover a weak-mode solution to 1e-12,
    # independently of RHS scale; permittivity-weighted bounds use exact sums.
    add_executable(test_gmres test_gmres.cpp)
    target_link_libraries(test_gmres FasterCapSolver)
    add_test(NAME gmres COMMAND test_gmres)
    set_tests_properties(gmres PROPERTIES TIMEOUT 10)
endif()
