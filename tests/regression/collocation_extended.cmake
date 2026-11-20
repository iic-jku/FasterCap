# Two interfaces (1280 panels) and uneven triangulation (1112 panels) retain
# the 2% analytic budget of the existing depth-2 sphere tests. Reciprocity has
# the same separate 0.5% budget, including conductor collocation error.
# These isolate discretization from grouping; the high-contrast hierarchy
# limitation remains documented and must not be hidden by these direct tests.
foreach(case layered_spheres_flux coated_spheres_flux_uneven)
    fastercap_test(collocation_${case} ${case} 0.02 RECIPROCITY_TOLERANCE 0.005
        -m1e32 -d1e-40 -t1e-9 -pj -f0)
    set_tests_properties(collocation_${case} PROPERTIES TIMEOUT 30)
endforeach()
