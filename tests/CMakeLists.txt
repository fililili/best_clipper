CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    GIT_TAG v1.17.0
    EXCLUDE_FROM_ALL YES
)

if(TARGET GTest::gtest_main)
    message(STATUS "GTest::gtest_main target found")
else()
    message(FATAL_ERROR "GTest::gtest_main target not found!")
endif()

add_executable(test_known_wkt test_known_wkt.cpp)
target_compile_features(test_known_wkt PRIVATE cxx_std_20)
target_link_libraries(test_known_wkt PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_known_wkt COMMAND test_known_wkt)

add_executable(test_connected_components test_connected_components.cpp)
target_compile_features(test_connected_components PRIVATE cxx_std_20)
target_link_libraries(test_connected_components PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_connected_components COMMAND test_connected_components)

add_executable(test_manhattan test_manhattan.cpp)
target_compile_features(test_manhattan PRIVATE cxx_std_20)
target_link_libraries(test_manhattan PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_manhattan COMMAND test_manhattan)

add_executable(test_uniform_grid test_uniform_grid.cpp)
target_compile_features(test_uniform_grid PRIVATE cxx_std_20)
target_link_libraries(test_uniform_grid PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_uniform_grid COMMAND test_uniform_grid)

add_executable(test_union_rect test_union_rect.cpp)
target_compile_features(test_union_rect PRIVATE cxx_std_20)
target_link_libraries(test_union_rect PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_union_rect COMMAND test_union_rect)

set_tests_properties(
    _test_known_wkt _test_union_rect
    _test_connected_components _test_manhattan _test_uniform_grid
    PROPERTIES LABELS "unit")
