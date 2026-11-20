# #78: -a0.1 must keep reciprocity, closed-conductor shielding and agreement with
# a tightly solved identical mesh within its requested 10% accuracy. Planar
# surfaces avoid input-polyhedron versus analytic-sphere ambiguity. Galerkin
# isolates automatic solve accuracy from the separate collocation defect #74.
# A tight solve on matched panels/links additionally checks the 2.5% algebraic
# sub-budget, separately from the 10% physical/direct-mesh envelope.
# Nine serial solves, at most 2048 panels, with 20 seconds per solver process.
add_test(NAME automatic_accuracy
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_automatic_accuracy.py
        --fastercap $<TARGET_FILE:FasterCap>)
set_tests_properties(automatic_accuracy PROPERTIES TIMEOUT 60)
