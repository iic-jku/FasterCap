# minimum version is 3.0 to be able to handle Eigen library
cmake_minimum_required (VERSION 3.12)
project (FasterCap)

# main project specific settings
#
# option to compile as headless, or with GUI
set(FASTFIELDSOLVERS_HEADLESS OFF CACHE BOOL "Compile as headless (i.e. no GUI)")

#
# Global Compiler and Linker Settings
#

#
# using only GCC, therefore configuring build only for GCC
#
# .cmake files are script files; not using here any script
#include(cMake/custom_include.cmake)
#configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
#
# if headless, pass on the "#define FCG_HEADLESS"
if(${FASTFIELDSOLVERS_HEADLESS} MATCHES ON)
    add_definitions("-DFCG_HEADLESS -DwxUSE_GUI=0")
else()
    add_definitions("-DwxUSE_GUI=1")
endif()
#
# pass to gcc the argument -DHAVE_CONFIG_H . This defines the pre-processor token HAVE_CONFIG_H
# exactly as if you had #define HAVE_CONFIG_H right at the start of each of your source files.
# HAVE_CONFIG_H is a pre-defined macro in the build-system that generates a config.h
# so that code knows whether it needs to #include config.h
#add_definitions(-DHAVE_CONFIG_H)
#

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set (EXTRA_LIBS ${EXTRA_LIBS} OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "Unable to find OpenMP. Try exporting the OpenMP_ROOT environmental variable. E.g., on macOS: export OpenMP_ROOT=$(brew --prefix)/opt/libomp")
endif()

# BLAS (library, cblas.h, and on macOS Accelerate) comes with the LinAlgebra target

if (NOT MSVC)
  set(FASTFIELDSOLVERS_COMPILE_OPTIONS ${FASTFIELDSOLVERS_COMPILE_OPTIONS} -Wextra -Wall)
endif()

# no fused multiply-add contraction, also in LinAlgebra and Geometry: Clang contracts by default
# on arm64, while x86-64 has no FMA instructions without -march, and -a can then stop after a
# different number of refinements (iic-jku/FasterCap#10)
if (NOT MSVC)
  add_compile_options(-ffp-contract=off)
endif()

# -s: strip all symbols from the executable
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")

# wxWidgets usage
# See CMake docs "FindwxWidgets"
# In general, for *nix CMake already generates the 'wx-confg' options that are needed
# for compiling and linking with wxWidgets. However, if anything special is needed,
# it must be tweaked into "wxWidgets_CONFIG_OPTIONS".
# In this case, we pass the options `wx-config --version=3.0 --static=no --debug`,
# but '--debug' is only needed in the debug build. Therefore, we use a CMake generator-expression
# to define the flag "--debug" only in case of debug build (this is nested: the first clause
# $<CONFIG:Debug> generates 1 or 0 if CONFIG is Debug or not. Then we use the conditional
# expression $<condition:true_string>, so if the previous evaluates to 1, here the string is "--debug",
# otherwise nothing.
set(wxWidgets_CONFIG_OPTIONS --static=no $<$<CONFIG:Debug>:--debug>)
# these flags are BUILD_TYPE dependent
# -D__WXDEBUG__: define that wxWidgets are in debug version
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D__WXDEBUG__")
if(${FASTFIELDSOLVERS_HEADLESS} MATCHES ON)
    # we only need the 'base' wxWidget package for headless compile
    find_package(wxWidgets COMPONENTS base)
else()
    # probably not all required. Note however that order is important.
    #find_package(wxWidgets COMPONENTS html xrc net core base msw png jpeg tiff zlib kernel32 user32 gdi32 winspool comdlg32 advapi32 shell32 ole32 oleaut32 uuid comctl32 sock32 odbc32 iphlpapi)
    find_package(wxWidgets COMPONENTS core base html)
endif()
# Configuration to use
# (e.g., msw, mswd, mswu, mswunivud, etc.)
# wxWidgets_CONFIGURATION
# general include (needed for wxWidgets in CMake), see CMake docs "FindwxWidgets"
include(${wxWidgets_USE_FILE})
# add to WXWIDGETS_LIBS the wxWidgets libraries
# we'll subsequently link WXWIDGETS_LIBS with 'target_link_libraries'
set(WXWIDGETS_LIBS ${WXWIDGETS_LIBS} ${wxWidgets_LIBRARIES})

# general include_directories, as LinAlgebra and Geometry
# are one level under the base path, and they are included
# in the sources as "#include "LinAlgebra/Mtx.h" etc.
# '..' resolved by CMake, not by the file system: with a symlinked FasterCap directory
# (e.g. a workspace of symlinks), the file system would resolve it to the symlink target's parent
get_filename_component(FASTFIELDSOLVERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
include_directories("${FASTFIELDSOLVERS_DIR}")
# same as above for the main includes in the base path; could use ${CMAKE_SOURCE_DIR}
include_directories(".")

# regression tests, run with: ctest --test-dir <build directory> --output-on-failure;
# enabled before adding LinAlgebra and Geometry, so that their tests are run too
enable_testing()

include_directories ("${FASTFIELDSOLVERS_DIR}/LinAlgebra")
add_subdirectory("${FASTFIELDSOLVERS_DIR}/LinAlgebra" LinAlgebra)
include_directories ("${FASTFIELDSOLVERS_DIR}/Geometry")
add_subdirectory("${FASTFIELDSOLVERS_DIR}/Geometry" Geometry)
set (EXTRA_LIBS ${EXTRA_LIBS} LinAlgebra Geometry)

include_directories ("${PROJECT_SOURCE_DIR}/Solver")
include("./Solver/CMakeLists.txt")
# do not create a library, to speed up compilation in case of file modifications
# (only obj to recompile and final linking, no need to re-create also the library)
#add_subdirectory(Solver)
#set (EXTRA_LIBS ${EXTRA_LIBS} Solver)

# the solver, compiled once for FasterCap and for the unit tests in tests/
add_library(FasterCapSolver OBJECT ${FasterCap_Sources} FasterCapGlobal.cpp FasterCapGlobal.h)
target_compile_options(FasterCapSolver PUBLIC ${FASTFIELDSOLVERS_COMPILE_OPTIONS})
target_link_libraries(FasterCapSolver PUBLIC ${EXTRA_LIBS} ${WXWIDGETS_LIBS} ${EIGEN_LIBS})

# The optional JSON report identifies the configuration and source state actually built.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
add_custom_target(FasterCapBuildProvenance
    COMMAND ${CMAKE_COMMAND}
        "-DFASTERCAP_SOURCE=${PROJECT_SOURCE_DIR}"
        "-DLINALGEBRA_SOURCE=${FASTFIELDSOLVERS_DIR}/LinAlgebra"
        "-DGEOMETRY_SOURCE=${FASTFIELDSOLVERS_DIR}/Geometry"
        "-DTEMPLATE=${PROJECT_SOURCE_DIR}/cmake/BuildProvenance.h.in"
        "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/generated/BuildProvenance.h"
        -P "${PROJECT_SOURCE_DIR}/cmake/BuildProvenance.cmake"
    BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/generated/BuildProvenance.h"
    VERBATIM)
add_dependencies(FasterCapSolver FasterCapBuildProvenance)
get_target_property(FC_REPORT_BLAS LinAlgebra INTERFACE_LINK_LIBRARIES)
if(NOT FC_REPORT_BLAS)
    set(FC_REPORT_BLAS "")
endif()
set(FC_REPORT_COMPILER "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
set(FC_REPORT_SYSTEM "${CMAKE_SYSTEM_NAME}")
set(FC_REPORT_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}")
# C++ string literals, including Windows library paths.
foreach(value FC_REPORT_COMPILER FC_REPORT_SYSTEM FC_REPORT_PROCESSOR FC_REPORT_BLAS)
    string(REPLACE "\\" "\\\\" ${value} "${${value}}")
    string(REPLACE "\"" "\\\"" ${value} "${${value}}")
endforeach()
configure_file(cmake/BuildConfiguration.h.in generated/BuildConfiguration.h.in @ONLY)
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/$<CONFIG>/BuildConfiguration.h"
    INPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/BuildConfiguration.h.in")
target_include_directories(FasterCapSolver PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/generated" "${CMAKE_CURRENT_BINARY_DIR}/generated/$<CONFIG>")

set(FastFieldSolvers_Sources ${FastFieldSolvers_Sources}
    FasterCapConsole.cpp
    FasterCapConsole.h
    test.cpp
    test.h)

# if running with a GUI
if(${FASTFIELDSOLVERS_HEADLESS} MATCHES OFF)
    set(FastFieldSolvers_Sources ${FastFieldSolvers_Sources}
        AboutBox.cpp
        AboutBox.h
        FasterCapApp.cpp
        FasterCapApp.h
        FasterCapMain.cpp
        FasterCapMain.h
        LicenseBox.cpp
        LicenseBox.h
        RunDialog.cpp
        RunDialog.h
        resource.rc)
    # and if also running on Windows
    if(CMAKE_SYSTEM_NAME STREQUAL Windows)
        set(FastFieldSolvers_Sources ${FastFieldSolvers_Sources}
            AutomationHelper.cpp
            AutomationHelper.h)
    endif()
endif()

add_executable (FasterCap ${FastFieldSolvers_Sources})
target_link_libraries(FasterCap FasterCapSolver)

add_subdirectory(tests)
