# #76: preconditioning changes convergence, not the solved matrix. With 1e-9
# GMRES tolerance, allow 1e-6 entrywise disagreement for conditioning/roundoff.
# Existing 2% depth-2 sphere and 1% coax analytic budgets are retained.
# The fixed-mesh coax uses 512 segments per surface: direct discretization error
# is 0.93%, versus 1.83% for 256 without the existing automatic test's refinement.
# Galerkin sphere cases isolate block correctness from the separate default
# collocation discretization issue #74; coax exercises the 2D segment kernels.
# Nine serial solves, at most 1536 panels each, without automatic refinement.
add_test(NAME block_preconditioner_matrix
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_block_preconditioner.py
        --fastercap $<TARGET_FILE:FasterCap>)
set_tests_properties(block_preconditioner_matrix PROPERTIES TIMEOUT 60)
