# Regression tests: check_capacitance.py runs FasterCap on a geometry of geometries.py
# and checks the matrix (exit code, signs, symmetry, analytic values, symmetric entries).
# Tolerances are relative; the comment of each test says why its tolerance is enough.

# macOS: LinAlgebra uses the Accelerate headers, so FasterCap must link Accelerate's BLAS,
# not one found elsewhere (MacPorts OpenBLAS)
if(APPLE)
    add_test(NAME macos_blas_is_accelerate
        COMMAND ${CMAKE_COMMAND} -DEXECUTABLE=$<TARGET_FILE:FasterCap>
            -P ${CMAKE_CURRENT_SOURCE_DIR}/check_blas_linkage.cmake)
endif()

# -bv prints only the version; headless builds started a run without input file instead
add_test(NAME version COMMAND FasterCap -bv)
set_tests_properties(version PROPERTIES
    PASS_REGULAR_EXPRESSION "^FasterCap version [0-9]+\\.[0-9]+\\.[0-9]+"
    FAIL_REGULAR_EXPRESSION "Running|Error")

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "Python 3 not found: no regression tests")
    return()
endif()

# fastercap_test(<name> <geometry> <tolerance> [OUT_OF_CORE] [INTERRUPT] [NO_MATRIX]
#                [EXIT_CODE <code>] <FasterCap options>...)
function(fastercap_test name geometry tolerance)
    cmake_parse_arguments(PARSE_ARGV 3 arg "OUT_OF_CORE;INTERRUPT;NO_MATRIX" "EXIT_CODE" "")
    set(flags)
    if(arg_OUT_OF_CORE)
        list(APPEND flags --out-of-core)
    endif()
    if(arg_INTERRUPT)
        list(APPEND flags --interrupt)
    endif()
    if(arg_NO_MATRIX)
        list(APPEND flags --no-matrix)
    endif()
    if(DEFINED arg_EXIT_CODE)
        list(APPEND flags --exit-code ${arg_EXIT_CODE})
    endif()
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_capacitance.py
            --fastercap $<TARGET_FILE:FasterCap>
            --geometry ${geometry}
            --tolerance ${tolerance}
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/${name}
            ${flags}
            -- ${arg_UNPARSED_ARGUMENTS})
    # a hang fails the test instead of blocking the CI job
    set_tests_properties(${name} PROPERTIES TIMEOUT 300)
endfunction()

# default settings, no -a: 18 panels, a few milliseconds; C12 = C13 within 0.1 %, symmetric
# within 0.6 %
fastercap_test(three_cubes_default three_cubes 0.01)
# 0.2 % below the analytic value at -a0.001
fastercap_test(cube cube 0.01 -a0.001)
# the inscribed polyhedron of 1280 triangles is 0.3 % below the sphere
fastercap_test(sphere sphere 0.01 -a0.01)
# C12 = C13 and C22 = C33 within 0.06 %, the matrix symmetric within 0.06 %
fastercap_test(three_cubes three_cubes 0.01 -a0.01)
# Galerkin: C12 and C13 differ by 0.8 %, also at -t1e-4
fastercap_test(three_cubes_galerkin three_cubes 0.02 -a0.01 -g)
# 2.4M links (over two link chunks) at -m0.005 -d0.5; C12 = C13 within 0.02 %, symmetric
# within 0.01 %; the out-of-core matrix must equal the in-core one
fastercap_test(three_cubes_out_of_core three_cubes 0.01 OUT_OF_CORE -m0.005 -d0.5)
# block preconditioner with blocks of 2: blocks rooted at a leaf were never completed, and the
# run hung after writing past the block element array; C12 = C13 within 0.07 %, symmetric
# within 0.3 %
fastercap_test(three_cubes_block_precond three_cubes 0.01 -pb2)
# dielectric shell: links between dielectric panels and super panels next to them made it 16 %
# high. The collocation on the dielectric interface converges slowly: with every link direct (and
# in FastCap2) the input panels give 3.3 % above the analytic value, -a0.01 gives 2.5 %
fastercap_test(coated_sphere coated_sphere 0.05 -a0.01)
# -t1e-30 is out of reach of GMRES: exit code 75 (FC_CONVERGENCE_FAILED) with the matrix of the
# last iteration, as accurate as three_cubes_default
fastercap_test(three_cubes_not_converged three_cubes 0.01 EXIT_CODE 75 -t1e-30)
# a triangle with its vertices on a line has no self-potential, the matrix would be NaN: exit
# code 65 (FC_INVALID_RESULT)
fastercap_test(degenerate_triangle degenerate_triangle 0 NO_MATRIX EXIT_CODE 65)
# FasterCap can't read PATRAN neutral files: exit code 74 (FC_FILE_ERROR), not the matrix
# without the dielectric
fastercap_test(patran_dielectric patran_dielectric 0 NO_MATRIX EXIT_CODE 74)
# SIGINT while out-of-core: exit code 125 (FC_USER_BREAK), temporary files deleted; Windows has
# no SIGINT for a single process
if(NOT WIN32)
    fastercap_test(three_cubes_interrupted three_cubes 0 INTERRUPT -m0.005 -d0.5)
endif()
