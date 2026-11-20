if(FASTFIELDSOLVERS_HEADLESS)
    add_executable(count_benchmark_links ../tools/count_benchmark_links.cpp)
    target_link_libraries(count_benchmark_links FasterCapSolver)
endif()
