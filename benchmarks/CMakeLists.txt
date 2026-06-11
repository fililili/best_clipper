set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
CPMAddPackage(
    NAME benchmark
    GITHUB_REPOSITORY google/benchmark
    GIT_TAG v1.9.5
    EXCLUDE_FROM_ALL YES
)
if(TARGET benchmark::benchmark)
    message(STATUS "benchmark::benchmark target found")
else()
    message(FATAL_ERROR "benchmark::benchmark target not found!")
endif()

# Clipper2: disable unwanted sub-projects before adding
set(CLIPPER2_TESTS OFF CACHE BOOL "" FORCE)
set(CLIPPER2_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CLIPPER2_UTILS OFF CACHE BOOL "" FORCE)
CPMAddPackage(
    NAME Clipper2
    GITHUB_REPOSITORY AngusJohnson/Clipper2
    GIT_TAG main
    SOURCE_SUBDIR CPP
)
if(TARGET Clipper2)
    message(STATUS "Clipper2 target found")
else()
    message(FATAL_ERROR "Clipper2 target not found!")
endif()

add_executable(bench_ops bench_ops.cpp)
target_compile_features(bench_ops PRIVATE cxx_std_20)
target_link_libraries(bench_ops PRIVATE best_clipper_core benchmark::benchmark Clipper2)
add_test(NAME _bench_ops COMMAND bench_ops)
set_tests_properties(_bench_ops PROPERTIES LABELS "benchmark")

add_executable(bench_ops_comparison bench_ops_comparison.cpp)
target_compile_features(bench_ops_comparison PRIVATE cxx_std_20)
target_link_libraries(bench_ops_comparison PRIVATE best_clipper_core benchmark::benchmark Clipper2)
add_test(NAME _bench_ops_comparison COMMAND bench_ops_comparison)
set_tests_properties(_bench_ops_comparison PROPERTIES LABELS "benchmark")