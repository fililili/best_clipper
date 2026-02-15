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

add_executable(basic_test1 test1.cpp)
target_compile_features(basic_test1 PRIVATE cxx_std_20)
target_link_libraries(basic_test1 PRIVATE best_clipper_basic GTest::gtest_main)
add_test(NAME _basic_test1 COMMAND basic_test1)