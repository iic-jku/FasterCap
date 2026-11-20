# Pure report/decision tests. The actual four-platform plate captures run in CI
# and compare each pair under the budgets documented in platform-comparison.md.
add_test(NAME platform_comparison
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_platform_comparison.py)
set_tests_properties(platform_comparison PROPERTIES TIMEOUT 10)
