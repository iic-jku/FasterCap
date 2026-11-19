# Fixed, all-direct Galerkin meshes isolate dielectric field integration. Depth-2
# inscribed spheres have roughly 1% radius error, allowing 2% physical matrix error.
# Depth 3 reduces that geometric error to about 0.3%, within a 0.5% envelope.
# Reciprocity has a separate 0.5% budget for source/far quadrature, tightened to
# 0.1% on the finer mesh. The old field calculation gives 2-33% defects here.
# Offset surfaces have no spherical analytic formula but must remain reciprocal.
# No automatic mesh growth: 960/3840 panels, at most 14,741,760 directed links.
foreach(case coated_spheres_flux coated_spheres_flux_lossy coated_spheres_flux_high_contrast
             coated_spheres_flux_offset coated_spheres_flux_thin)
    fastercap_test(${case} ${case} 0.02 RECIPROCITY_TOLERANCE 0.005
        -g -m1e32 -d1e-40 -t1e-9 -pj -f0)
    set_tests_properties(${case} PROPERTIES TIMEOUT 30)
endforeach()
fastercap_test(coated_spheres_flux_fine coated_spheres_flux_fine 0.005
    RECIPROCITY_TOLERANCE 0.001 -g -m1e32 -d1e-40 -t1e-9 -pj -f0)
set_tests_properties(coated_spheres_flux_fine PROPERTIES TIMEOUT 30)
