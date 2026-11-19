if(UNIX)
    add_test(NAME resource_guard
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_resource_guard.py)
    set_tests_properties(resource_guard PROPERTIES TIMEOUT 30 RUN_SERIAL TRUE)
endif()
