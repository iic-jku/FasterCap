# Closed conductor fixtures need consistently outward face normals.
add_test(NAME geometry_fixtures COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_geometries.py)
