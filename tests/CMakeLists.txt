CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    GIT_TAG v1.17.0
)

if(TARGET GTest::gtest_main)
    message(STATUS "GTest::gtest_main target found")
else()
    message(FATAL_ERROR "GTest::gtest_main target not found!")
endif()

add_executable(test_union test_union.cpp)
target_compile_features(test_union PRIVATE cxx_std_20)
target_link_libraries(test_union PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_union COMMAND test_union)

add_executable(test_union_rect test_union_rect.cpp)
target_compile_features(test_union_rect PRIVATE cxx_std_20)
target_link_libraries(test_union_rect PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_union_rect COMMAND test_union_rect)

add_executable(test_intersection test_intersection.cpp)
target_compile_features(test_intersection PRIVATE cxx_std_20)
target_link_libraries(test_intersection PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_intersection COMMAND test_intersection)

add_executable(test_xor test_xor.cpp)
target_compile_features(test_xor PRIVATE cxx_std_20)
target_link_libraries(test_xor PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_xor COMMAND test_xor)

add_executable(test_difference test_difference.cpp)
target_compile_features(test_difference PRIVATE cxx_std_20)
target_link_libraries(test_difference PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_difference COMMAND test_difference)

add_executable(test_connected_components test_connected_components.cpp)
target_compile_features(test_connected_components PRIVATE cxx_std_20)
target_link_libraries(test_connected_components PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_connected_components COMMAND test_connected_components)

add_executable(test_manhattan test_manhattan.cpp)
target_compile_features(test_manhattan PRIVATE cxx_std_20)
target_link_libraries(test_manhattan PRIVATE best_clipper_core GTest::gtest_main)
add_test(NAME _test_manhattan COMMAND test_manhattan)
