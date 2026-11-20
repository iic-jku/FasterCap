# Coarse inscribed spheres retain the existing 2% analytic envelope. Allow only
# 0.5% additional matrix error from grouping at fixed mesh size, and 1% total
# reciprocity error including the roughly 0.4% conductor-collocation defect.
# This is a moderate-contrast regression, not a high-contrast accuracy guarantee.
# Nine serial solves, exactly 960 panels, no automatic mesh growth.
add_test(NAME dielectric_hierarchy_matrix
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_dielectric_hierarchy.py
        --fastercap $<TARGET_FILE:FasterCap>
        --workdir ${CMAKE_CURRENT_BINARY_DIR}/dielectric_hierarchy_matrix)
set_tests_properties(dielectric_hierarchy_matrix PROPERTIES TIMEOUT 60)

# Unequal panel areas and multiple interfaces retain the same 2% analytic and
# 1% reciprocity budgets. Limit additional compression error to 0.2% at this
# tight interaction setting, one tenth of the mesh's physical-error envelope.
# Six serial solves, at most 1280 panels, including direct fixed-mesh oracles.
add_test(NAME dielectric_hierarchy_uneven
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_dielectric_hierarchy.py
        --fastercap $<TARGET_FILE:FasterCap> --uneven
        --workdir ${CMAKE_CURRENT_BINARY_DIR}/dielectric_hierarchy_uneven)
set_tests_properties(dielectric_hierarchy_uneven PROPERTIES TIMEOUT 60)

# At contrast 100, hierarchy errors are amplified by the dielectric solve.
# Budget 1.5% additional error relative to the identical direct mesh, while
# retaining the 2% analytic and 1% reciprocity limits. The coarse direct mesh
# itself has about 1% geometric error. Nine serial 960-panel solves exercise
# concentric, eccentric and thin layers without automatic mesh growth.
add_test(NAME dielectric_hierarchy_high_contrast
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_dielectric_hierarchy.py
        --fastercap $<TARGET_FILE:FasterCap>
        --workdir ${CMAKE_CURRENT_BINARY_DIR}/dielectric_hierarchy_high_contrast
        --high-contrast)
set_tests_properties(dielectric_hierarchy_high_contrast PROPERTIES TIMEOUT 60)
