# Opt-in JSON: exact double round trips, complex/2D values, valid failure reports,
# and identical legacy console/CSV matrices with and without the report.
if(FASTFIELDSOLVERS_HEADLESS)
    add_executable(test_numerical_report test_numerical_report.cpp)
    target_link_libraries(test_numerical_report FasterCapSolver)
    add_test(NAME numerical_report_json
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_numerical_report.py
            --writer $<TARGET_FILE:test_numerical_report>
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/numerical_report_json)
endif()
add_test(NAME numerical_report_solver
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_numerical_report.py
        --fastercap $<TARGET_FILE:FasterCap>
        --workdir ${CMAKE_CURRENT_BINARY_DIR}/numerical_report_solver)
set_tests_properties(numerical_report_solver PROPERTIES TIMEOUT 180)
