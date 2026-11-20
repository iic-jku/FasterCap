# The default collocation scheme uses the same physical bounds as the Galerkin probe.
# Depth-2 geometric errors are about 1%; depth 3 about 0.3%. Maximum 3840 panels.
foreach(case coated_spheres_flux coated_spheres_flux_lossy coated_spheres_flux_high_contrast
             coated_spheres_flux_offset coated_spheres_flux_thin)
    fastercap_test(collocation_${case} ${case} 0.02 RECIPROCITY_TOLERANCE 0.005
        -m1e32 -d1e-40 -t1e-9 -pj -f0)
    set_tests_properties(collocation_${case} PROPERTIES TIMEOUT 30)
endforeach()
fastercap_test(collocation_coated_spheres_flux_fine coated_spheres_flux_fine 0.005
    RECIPROCITY_TOLERANCE 0.001 -m1e32 -d1e-40 -t1e-9 -pj -f0)
set_tests_properties(collocation_coated_spheres_flux_fine PROPERTIES TIMEOUT 30)
