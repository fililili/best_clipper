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

add_executable(bench_comparison_rect bench_comparison_rect.cpp)
target_compile_features(bench_comparison_rect PRIVATE cxx_std_20)
target_link_libraries(bench_comparison_rect PRIVATE best_clipper_core benchmark::benchmark Clipper2)
add_test(NAME _bench_comparison_rect COMMAND bench_comparison_rect)
set_tests_properties(_bench_comparison_rect PROPERTIES LABELS "benchmark")

add_executable(bench_comparison_circle bench_comparison_circle.cpp)
target_compile_features(bench_comparison_circle PRIVATE cxx_std_20)
target_link_libraries(bench_comparison_circle PRIVATE best_clipper_core benchmark::benchmark Clipper2)
add_test(NAME _bench_comparison_circle COMMAND bench_comparison_circle)
set_tests_properties(_bench_comparison_circle PROPERTIES LABELS "benchmark")