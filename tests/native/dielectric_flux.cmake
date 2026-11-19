if(FASTFIELDSOLVERS_HEADLESS)
    # Gauss flux on a closed polyhedron is 1, 0 or 1/2 for a unit source inside,
    # outside or on a smooth face. Dyadic coordinates leave only roundoff (2e-12).
    add_executable(test_dielectric_flux test_dielectric_flux.cpp)
    target_link_libraries(test_dielectric_flux FasterCapSolver)
    add_test(NAME dielectric_flux COMMAND test_dielectric_flux)
    set_tests_properties(dielectric_flux PROPERTIES TIMEOUT 10)
endif()
