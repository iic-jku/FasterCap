# #85: required C/D fields must fail at their source line, before recursion or
# numerical solving. Valid spelling/marker variants must dump identical geometry.
add_test(NAME material_records
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_material_records.py
        $<TARGET_FILE:FasterCap>)
set_tests_properties(material_records PROPERTIES TIMEOUT 60)
