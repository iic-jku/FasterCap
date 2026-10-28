# Regression tests: check_capacitance.py runs FasterCap on a geometry of geometries.py
# and checks the matrix (exit code, signs, symmetry, analytic values, symmetric entries).
# Tolerances are relative; the comment of each test says why its tolerance is enough.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "Python 3 not found: no regression tests")
    return()
endif()

function(fastercap_test name geometry tolerance)
    add_test(NAME ${name}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_capacitance.py
            --fastercap $<TARGET_FILE:FasterCap>
            --geometry ${geometry}
            --tolerance ${tolerance}
            --workdir ${CMAKE_CURRENT_BINARY_DIR}/${name}
            -- ${ARGN})
endfunction()

# 0.2 % below the analytic value at -a0.001
fastercap_test(cube cube 0.01 -a0.001)
# the inscribed polyhedron of 1280 triangles is 0.3 % below the sphere
fastercap_test(sphere sphere 0.01 -a0.01)
# C12 = C13 and C22 = C33 within 0.06 %, the matrix symmetric within 0.06 %
fastercap_test(three_cubes three_cubes 0.01 -a0.01)
# Galerkin: C12 and C13 differ by 0.8 %, also at -t1e-4
fastercap_test(three_cubes_galerkin three_cubes 0.02 -a0.01 -g)
