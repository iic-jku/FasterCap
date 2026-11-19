# Exercise the production discovery helper without compiling FasterCap or needing Python.
file(REMOVE_RECURSE "${TEST_ROOT}/source" "${TEST_ROOT}/build")
set(source "${TEST_ROOT}/source")
set(build "${TEST_ROOT}/build")
file(MAKE_DIRECTORY "${source}/registrations" "${build}")
file(WRITE "${source}/CMakeLists.txt"
    "cmake_minimum_required(VERSION 3.12)\n"
    "project(TestRegistration NONE)\n"
    "enable_testing()\n"
    "include(\"${REGISTRATION_HELPER}\")\n"
    "fastercap_include_tests(registrations)\n")
# Scratch files must not become registration modules.
file(WRITE "${source}/registrations/notes.txt" "not CMake code\n")

function(run_checked)
    execute_process(COMMAND ${ARGV} RESULT_VARIABLE result
        OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 20
        WORKING_DIRECTORY "${build}")
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "${ARGV} failed (${result}):\n${output}\n${error}")
    endif()
endfunction()

set(generator_args -G "${TEST_GENERATOR}")
if(TEST_GENERATOR_PLATFORM)
    list(APPEND generator_args -A "${TEST_GENERATOR_PLATFORM}")
endif()
if(TEST_GENERATOR_TOOLSET)
    list(APPEND generator_args -T "${TEST_GENERATOR_TOOLSET}")
endif()
run_checked("${CMAKE_COMMAND}" "${source}" ${generator_args})

function(expect_tests)
    # Only build after changing the module set: no explicit reconfigure is allowed here.
    run_checked("${CMAKE_COMMAND}" --build "${build}" --config "${TEST_CONFIG}" --parallel 1)
    execute_process(COMMAND "${CMAKE_CTEST_COMMAND}" -N -C "${TEST_CONFIG}"
        WORKING_DIRECTORY "${build}" RESULT_VARIABLE result
        OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 20)
    if(NOT "${result}" STREQUAL "0")
        message(FATAL_ERROR "CTest discovery failed (${result}):\n${output}\n${error}")
    endif()
    string(REGEX MATCHALL "Test[ \t]+#[0-9]+: [^\r\n]+" entries "${output}")
    set(names)
    foreach(entry IN LISTS entries)
        string(REGEX REPLACE "^Test[ \t]+#[0-9]+: " "" name "${entry}")
        list(APPEND names "${name}")
    endforeach()
    if(NOT "${names}" STREQUAL "${ARGV}")
        message(FATAL_ERROR "Expected tests '${ARGV}', got '${names}':\n${output}")
    endif()
endfunction()

expect_tests()
# Create in reverse order; discovery must be deterministic and include both additions.
file(WRITE "${source}/registrations/zeta.cmake"
    "add_test(NAME zeta COMMAND \"${CMAKE_COMMAND}\" -E true)\n")
file(WRITE "${source}/registrations/alpha.cmake"
    "add_test(NAME alpha COMMAND \"${CMAKE_COMMAND}\" -E true)\n")
expect_tests(alpha zeta)
file(REMOVE "${source}/registrations/alpha.cmake")
expect_tests(zeta)
file(REMOVE "${source}/registrations/zeta.cmake")
expect_tests()
