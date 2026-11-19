# Construction/destruction and cleanup before any input need no disk blocks.
if(FASTFIELDSOLVERS_HEADLESS)
    add_executable(test_refiner_initialization test_refiner_initialization.cpp)
    target_link_libraries(test_refiner_initialization FasterCapSolver)
    add_test(NAME refiner_initialization COMMAND test_refiner_initialization)
endif()
