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

# unit tests of the solver code; in GUI builds the solver's messages need the GUI application
if(FASTFIELDSOLVERS_HEADLESS)
    add_executable(test_potential test_potential.cpp)
    target_link_libraries(test_potential FasterCapSolver)
    add_test(NAME potential COMMAND test_potential)
    add_executable(test_triangulate test_triangulate.cpp)
    target_link_libraries(test_triangulate FasterCapSolver)
    add_test(NAME triangulate COMMAND test_triangulate)
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
#                [DUMP_INPUT] [OUTPUT_FILES] [EXIT_CODE <code>] [REFERENCE <geometry>]
#                [NAMES <name,...>] [SAME_AS <other options>] <FasterCap options>...)
function(fastercap_test name geometry tolerance)
    cmake_parse_arguments(PARSE_ARGV 3 arg "OUT_OF_CORE;INTERRUPT;NO_MATRIX;DUMP_INPUT;OUTPUT_FILES"
                          "EXIT_CODE;REFERENCE;NAMES;SAME_AS" "")
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
    if(arg_DUMP_INPUT)
        list(APPEND flags --dump-input)
    endif()
    if(arg_OUTPUT_FILES)
        list(APPEND flags --output-files)
    endif()
    if(DEFINED arg_EXIT_CODE)
        list(APPEND flags --exit-code ${arg_EXIT_CODE})
    endif()
    if(DEFINED arg_REFERENCE)
        list(APPEND flags --reference ${arg_REFERENCE})
    endif()
    if(DEFINED arg_NAMES)
        list(APPEND flags --names ${arg_NAMES})
    endif()
    if(DEFINED arg_SAME_AS)
        list(APPEND flags "--same-as=${arg_SAME_AS}")
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
# #7: both inscribed polyhedra are about 0.3 % below their spheres, so is the matrix; symmetric
# within 0.1 %, and row 1 sums to zero (no capacitance to infinity) within 0.1 % of C11
fastercap_test(concentric_spheres concentric_spheres 0.01 -a0.001)
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
# a non-planar Q panel is split into the two triangles along its smaller diagonal; without -v,
# it went to the planar triangulator, 6 % lower
fastercap_test(skewed_quad skewed_quad 0.01 REFERENCE skewed_quad_triangles)
# a self-intersecting Q panel: exit code 74 (FC_FILE_ERROR) naming its line; its triangulation
# used freed memory and crashed
fastercap_test(bowtie_quad bowtie_quad 0 NO_MATRIX EXIT_CODE 74)
# 2D coaxial cable with a dielectric coat (er 10): links between dielectric segments and super
# segments next to them made it 6.7 % high. The polygons of 256 segments are 0.01 % below the
# circles; -a0.01 gives +0.24 %, +1.1 % with the 3D distance factor of 2
fastercap_test(coated_coax coated_coax 0.01 -a0.01)
# Jacobi preconditioner (-pj): C12 = C13 and C22 = C33 within 0.1 %, symmetric within 0.1 %
fastercap_test(three_cubes_jacobi three_cubes 0.01 -pj -a0.01)
# two-levels preconditioner of dimension 16 (-ps16): C12 = C13 within 0.2 %
fastercap_test(three_cubes_two_levels three_cubes 0.01 -ps16 -a0.01)
# mesh curvature coefficient (-mc1): the sphere is 0.03 % below its analytic value, as the
# inscribed polyhedron is 0.3 % below the sphere
fastercap_test(sphere_curvature sphere 0.01 -a0.01 -mc1)
# '+' merges two conductors into one: the same panels as one conductor in one file
fastercap_test(two_cubes_merged two_cubes_merged 0.01 REFERENCE two_cubes -a0.01)
# the input geometry dumped by -oi gives the same matrix
fastercap_test(three_cubes_dump_input three_cubes 0.01 DUMP_INPUT -a0.01)
# -o, -e and -c write the refined panels, the printed matrix and charge densities that add up
# to the matrix
fastercap_test(three_cubes_output_files three_cubes 0.01 OUTPUT_FILES -a0.01)
# a single input file with the cube as an internal 'File' section: the same panels
fastercap_test(three_cubes_single_file three_cubes_single_file 0.01 REFERENCE three_cubes -a0.01)
# 'N' renames the conductor of the cube's panels
fastercap_test(three_cubes_renamed three_cubes_renamed 0.01 NAMES g1_box,g2_box,g3_box REFERENCE three_cubes -a0.01)
# 'N' to an existing name merges the two conductors
fastercap_test(two_cubes_named two_cubes_named 0.01 REFERENCE two_cubes -a0.01)
# reference points outside each panel of the coat, as klayout-pex writes them, instead of one
# inside: the same dielectric sides
fastercap_test(coated_sphere_local_reference coated_sphere_local_reference 0.05 REFERENCE coated_sphere -a0.01)
# lossy coat (er 4-j1): 1.6 % from the analytic complex value, of which 0.8 % in the imaginary
# part, as the lossless coat gives 2.5 % at -a0.01
fastercap_test(coated_sphere_lossy coated_sphere_lossy 0.05 -a0.01)
# -r, -i and -v only print more: the same matrix
fastercap_test(three_cubes_logging three_cubes 0.01 SAME_AS -a0.01 -a0.01 -r -i -v)
# options removed in #31 stop with exit code 64 (FC_COMMAND_LINE_ERROR) instead of being ignored
# (-s, -kc, -km) or hanging (-ph)
foreach(option -s -s0.1 -kc -km -ph -phs1 -phe0.1 -pht0.5)
    fastercap_test(removed_option${option} three_cubes 0 NO_MATRIX EXIT_CODE 64 ${option})
endforeach()
