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

add_executable(bench_ops bench_ops.cpp)
target_compile_features(bench_ops PRIVATE cxx_std_20)
target_link_libraries(bench_ops PRIVATE best_clipper_core benchmark::benchmark)
add_test(NAME _bench_ops COMMAND bench_ops)
set_tests_properties(_bench_ops PROPERTIES LABELS "benchmark")