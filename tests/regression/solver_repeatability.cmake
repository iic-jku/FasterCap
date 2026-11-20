# Identical work must give identical double matrices, including repeated one-
# thread processes with different allocation addresses. Fixed 960-panel meshes.
add_test(NAME solver_repeatability
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_solver_repeatability.py
        --fastercap $<TARGET_FILE:FasterCap>)
set_tests_properties(solver_repeatability PROPERTIES TIMEOUT 45)
