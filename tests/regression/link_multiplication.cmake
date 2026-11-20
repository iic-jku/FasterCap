# The fixed 1280-panel dielectric mesh spans two 2^20-link chunks. Exact raw
# matrix/work equality detects skipped, repeated or reordered block contributions;
# this test introduces no numerical tolerance or performance threshold.
add_test(NAME link_multiplication_storage
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_link_multiplication.py
        --fastercap $<TARGET_FILE:FasterCap>)
set_tests_properties(link_multiplication_storage PROPERTIES TIMEOUT 90)
