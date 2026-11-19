add_test(NAME test_registration
    COMMAND ${CMAKE_COMMAND}
        "-DREGISTRATION_HELPER=${CMAKE_CURRENT_SOURCE_DIR}/RegisterTests.cmake"
        "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/test_registration"
        "-DTEST_GENERATOR=${CMAKE_GENERATOR}"
        "-DTEST_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}"
        "-DTEST_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}"
        "-DTEST_CONFIG=$<CONFIG>"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_registration.cmake)
set_tests_properties(test_registration PROPERTIES TIMEOUT 60)
