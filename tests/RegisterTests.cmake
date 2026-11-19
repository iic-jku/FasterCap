# Each concern owns its registration module; independent additions need no shared list.
# CONFIGURE_DEPENDS also notices modules added or removed after the initial configure.
# Include native modules before Python discovery, regression modules after its guard.
function(fastercap_include_tests directory)
    file(GLOB modules LIST_DIRECTORIES FALSE CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/${directory}/*.cmake")
    foreach(module IN LISTS modules)
        include("${module}")
    endforeach()
endfunction()
