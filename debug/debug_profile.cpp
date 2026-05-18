#define BEST_CLIPPER_PROFILE 1
#include <cstdio>
#include <cstdlib>
#include <random>
#include <chrono>
#include "core.hpp"
using namespace best_clipper;

ring make_rect(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2) {
    ring r;
    r.push_back(point{x1, y1});
    r.push_back(point{x1, y2});
    r.push_back(point{x2, y2});
    r.push_back(point{x2, y1});
    r.push_back(point{x1, y1});
    return r;
}

int main(int argc, char** argv) {
    size_t N = 1000000;
    if (argc > 1) N = (size_t)std::atoll(argv[1]);

    // Uniformly distributed rectangles across 200000 x 200000 plane
    // Rectangle sizes 2~30 in x and y, density ensures local overlap
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<uint32_t> xy(0, 200000);
    std::uniform_int_distribution<uint32_t> wh(2, 30);

    multi_polygon a;
    a.reserve(N);
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < N; i++) {
            uint32_t x = xy(rng), y = xy(rng);
            uint32_t w = wh(rng), h = wh(rng);
            a.push_back({{make_rect(x, y, x + w, y + h)}});
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        std::printf("=== Generated %zu rectangles in %.1f ms ===\n\n", N,
            std::chrono::duration<double, std::milli>(t1 - t0).count());
    }

    std::printf("=== Profiling self_or: %zu polygons ===\n", a.size());
    auto result = self_or(a);
    std::printf("Result: %zu polygons, area=%.0f\n\n", result.size(), bg::area(result));

    return 0;
}
