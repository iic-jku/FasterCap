# The same coarse inscribed spheres retain the existing 2% analytic and 1%
# reciprocity budgets. Allow only 0.2% additional error from conductor compression,
# one tenth of the mesh's physical-error envelope. The unbounded hierarchy adds
# 3.7% entry-relative error on this fixed mesh; physical bounds reduce it below 0.1%.
# Two serial solves, exactly 640 panels, at most 408960 directed physical links.
add_test(NAME conductor_hierarchy_matrix
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_dielectric_hierarchy.py
        --fastercap $<TARGET_FILE:FasterCap> --conductor
        --workdir ${CMAKE_CURRENT_BINARY_DIR}/conductor_hierarchy_matrix)
set_tests_properties(conductor_hierarchy_matrix PROPERTIES TIMEOUT 30)
