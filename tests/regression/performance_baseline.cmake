# Timing is opt-in and must preserve full-precision results. This checks manual
# and automatic solves, never a machine-dependent speed threshold.
if(TARGET count_benchmark_links)
    add_test(NAME performance_timings
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_performance_timings.py
            --fastercap $<TARGET_FILE:FasterCap>
            --counter $<TARGET_FILE:count_benchmark_links>)
    set_tests_properties(performance_timings PROPERTIES TIMEOUT 45)
endif()
