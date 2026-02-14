CPMAddPackage(
    NAME benchmark
    GITHUB_REPOSITORY google/benchmark
    GIT_TAG v1.9.5
)
if(TARGET benchmark::benchmark)
    message(STATUS "benchmark::benchmark target found")
else()
    message(FATAL_ERROR "benchmark::benchmark target not found!")
endif()

add_executable(basic_benchmark1 benchmark1.cpp)
target_compile_features(basic_benchmark1 PRIVATE cxx_std_20)
target_link_libraries(basic_benchmark1 PRIVATE best_clipper_basic)
add_test(NAME _basic_benchmark1 COMMAND basic_benchmark1)
set_tests_properties(_basic_benchmark1 PROPERTIES LABELS "benchmark")