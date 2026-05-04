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

add_executable(basic_test_union test_union.cpp)
target_compile_features(basic_test_union PRIVATE cxx_std_20)
target_link_libraries(basic_test_union PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_union COMMAND basic_test_union)

add_executable(basic_test_union_rect test_union_rect.cpp)
target_compile_features(basic_test_union_rect PRIVATE cxx_std_20)
target_link_libraries(basic_test_union_rect PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_union_rect COMMAND basic_test_union_rect)

add_executable(basic_test_intersection test_intersection.cpp)
target_compile_features(basic_test_intersection PRIVATE cxx_std_20)
target_link_libraries(basic_test_intersection PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_intersection COMMAND basic_test_intersection)

add_executable(basic_test_xor test_xor.cpp)
target_compile_features(basic_test_xor PRIVATE cxx_std_20)
target_link_libraries(basic_test_xor PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_xor COMMAND basic_test_xor)

add_executable(basic_test_difference test_difference.cpp)
target_compile_features(basic_test_difference PRIVATE cxx_std_20)
target_link_libraries(basic_test_difference PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_difference COMMAND basic_test_difference)

add_executable(basic_test_connected_components test_connected_components.cpp)
target_compile_features(basic_test_connected_components PRIVATE cxx_std_20)
target_link_libraries(basic_test_connected_components PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_connected_components COMMAND basic_test_connected_components)

add_executable(basic_test_manhattan test_manhattan.cpp)
target_compile_features(basic_test_manhattan PRIVATE cxx_std_20)
target_link_libraries(basic_test_manhattan PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test_manhattan COMMAND basic_test_manhattan)
