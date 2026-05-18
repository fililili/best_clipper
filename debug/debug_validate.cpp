// Quick validation: check best_clipper vs Boost area agreement
// using the same polygon generator as the benchmark.
#include "uint32_adaptor.hpp"
#include <random>
#include <cstdio>
#include <cmath>
#include <algorithm>

std::mt19937 rng;
constexpr int kSides = 11;

polygon_s32 gen_poly(int32_t cx, int32_t cy, int32_t radius, std::mt19937& rng) {
    std::uniform_int_distribution<int32_t> jitter(-2, 2);
    polygon_s32 poly;
    auto& outer = poly.outer();
    for (int j = 0; j < kSides; j++) {
        double a = 2.0 * 3.141592653589793 * j / kSides;
        int32_t vx = cx + (int32_t)(radius * std::cos(a)) + jitter(rng);
        int32_t vy = cy + (int32_t)(radius * std::sin(a)) + jitter(rng);
        outer.push_back({vx, vy});
    }
    std::reverse(outer.begin(), outer.end());
    outer.push_back(outer.front());
    poly.outer() = std::move(outer);
    return poly;
}

int main() {
    std::uniform_int_distribution<int32_t> radius_dist(8, 20);
    std::uniform_int_distribution<int32_t> jitter(-15, 15);

    for (int n : {10, 50, 100}) {
        // Same grid layout as RandPolygonsGrid in benchmark
        int cols = (int)std::ceil(std::sqrt((double)n));
        auto gen_set = [&](uint32_t seed) {
            rng.seed(seed);
            multi_polygon_s32 mp;
            for (int i = 0; i < n; i++) {
                int row = i / cols, col = i % cols;
                int32_t cx = 80 + col * 90 + jitter(rng);
                int32_t cy = 80 + row * 90 + jitter(rng);
                mp.push_back(gen_poly(cx, cy, radius_dist(rng), rng));
            }
            return mp;
        };
        auto a = gen_set(42);
        auto b = gen_set(99);

        int va = bg::is_valid(a), vb = bg::is_valid(b);
        printf("n=%d: a.mp_valid=%d b.mp_valid=%d\n", n, va, vb);

        multi_polygon_s32 bu, bi, bx, bd;
        bg::union_(a, b, bu);
        bg::intersection(a, b, bi);
        bg::sym_difference(a, b, bx);
        bg::difference(a, b, bd);
        auto cu = add(a, b);
        auto ci = intersection(a, b);
        auto cx = xor_(a, b);
        auto cd = difference(a, b);

        auto check = [](const char* op, double ba, double bc) {
            double diff = std::abs(ba - bc);
            double pct = 100.0 * diff / std::max(1.0, ba);
            printf("  %s: boost=%.0f bc=%.0f diff=%.0f (%.1f%%) %s\n",
                op, ba, bc, diff, pct, diff <= 1.0 ? "OK" : "");
        };
        check("union", bg::area(bu), bg::area(cu));
        check("inter", bg::area(bi), bg::area(ci));
        check("xor  ", bg::area(bx), bg::area(cx));
        check("diff ", bg::area(bd), bg::area(cd));
    }
    return 0;
}
