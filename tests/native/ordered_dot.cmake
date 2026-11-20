# These cancellation patterns have an exactly representable sum, independent of
# pointer alignment. Reject nonfinite terms rather than accepting a finite prefix.
add_executable(test_ordered_dot test_ordered_dot.cpp)
add_test(NAME ordered_dot COMMAND test_ordered_dot)
set_tests_properties(ordered_dot PROPERTIES TIMEOUT 10)
