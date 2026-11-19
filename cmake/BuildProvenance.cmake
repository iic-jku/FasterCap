# Run at build time so a checkout/edit after configuration does not leave stale source IDs.
find_package(Git QUIET)
foreach(repository FASTERCAP LINALGEBRA GEOMETRY)
    set(${repository}_REVISION "")
    set(${repository}_DIRTY null)
    if(GIT_FOUND AND EXISTS "${${repository}_SOURCE}/.git")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${${repository}_SOURCE}" rev-parse HEAD
            RESULT_VARIABLE result OUTPUT_VARIABLE revision OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(result EQUAL 0)
            set(${repository}_REVISION "${revision}")
            execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${${repository}_SOURCE}"
                status --porcelain --untracked-files=normal
                RESULT_VARIABLE result OUTPUT_VARIABLE changes ERROR_QUIET)
            if(result EQUAL 0)
                if(changes STREQUAL "")
                    set(${repository}_DIRTY false)
                else()
                    set(${repository}_DIRTY true)
                endif()
            endif()
        endif()
    endif()
endforeach()
configure_file("${TEMPLATE}" "${OUTPUT}" @ONLY)
