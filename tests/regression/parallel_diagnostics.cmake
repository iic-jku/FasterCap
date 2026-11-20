# This concerns diagnostics, not a different physical discretization. Keep
# existing geometry budgets and require exact full-precision matrix/count equality.
foreach(case three_cubes coated_coax near_triangles)
    add_test(NAME parallel_diagnostics_${case}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_parallel_diagnostics.py
            --fastercap $<TARGET_FILE:FasterCap> --geometry ${case})
    set_tests_properties(parallel_diagnostics_${case} PROPERTIES TIMEOUT 45)
endforeach()
