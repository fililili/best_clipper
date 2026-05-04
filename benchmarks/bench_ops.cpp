#include <benchmark/benchmark.h>
#include "uint32_adaptor.hpp"

// ============================================================================
// Fixture: pre-build N pairs of overlapping squares for union benchmarking
// ============================================================================
struct RectangleUnionFixture : benchmark::Fixture {
    void SetUp(const benchmark::State& state) override {
        int n = state.range(0);
        if (!a.empty() && (int)a.size() == n) return;
        a.clear(); b.clear();
        a.reserve(n); b.reserve(n);
        for (int i = 0; i < n; i++) {
            a.push_back(polygon_s32{{{{0+2*i,0+2*i},{0+2*i,2+2*i},{2+2*i,2+2*i},{2+2*i,0+2*i},{0+2*i,0+2*i}}}});
            b.push_back(polygon_s32{{{{1+2*i,1+2*i},{1+2*i,3+2*i},{3+2*i,3+2*i},{3+2*i,1+2*i},{1+2*i,1+2*i}}}});
        }
    }
    multi_polygon_s32 a, b;
};

// Small sizes: CI-friendly, multiple iterations
BENCHMARK_DEFINE_F(RectangleUnionFixture, RectangleUnion_Small)(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK_REGISTER_F(RectangleUnionFixture, RectangleUnion_Small)->Arg(10)->Arg(100)->Arg(1000);

// Large sizes: stress testing
BENCHMARK_DEFINE_F(RectangleUnionFixture, RectangleUnion_Large)(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK_REGISTER_F(RectangleUnionFixture, RectangleUnion_Large)
    ->Arg(10000)->Arg(100000)->Iterations(3);

// Million-scale: single-shot stress tests
BENCHMARK_DEFINE_F(RectangleUnionFixture, RectangleUnion_Million)(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK_REGISTER_F(RectangleUnionFixture, RectangleUnion_Million)
    ->Arg(500000)->Arg(1000000)->Iterations(1)->Unit(benchmark::kMillisecond);

// ============================================================================
// Fixture: pre-build N overlapping squares for self_union
// ============================================================================
struct RectangleSelfOrFixture : benchmark::Fixture {
    void SetUp(const benchmark::State& state) override {
        int n = state.range(0);
        if (!poly.empty() && (int)poly.size() == n) return;
        poly.clear();
        poly.reserve(n);
        for (int i = 0; i < n; i++)
            poly.push_back(polygon_s32{{{{0+i,0+i},{0+i,2+i},{2+i,2+i},{2+i,0+i},{0+i,0+i}}}});
    }
    multi_polygon_s32 poly;
};

BENCHMARK_DEFINE_F(RectangleSelfOrFixture, RectangleSelfOr_Small)(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(self_or(poly));
}
BENCHMARK_REGISTER_F(RectangleSelfOrFixture, RectangleSelfOr_Small)->Arg(10)->Arg(100)->Arg(1000);

BENCHMARK_DEFINE_F(RectangleSelfOrFixture, RectangleSelfOr_Large)(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(self_or(poly));
}
BENCHMARK_REGISTER_F(RectangleSelfOrFixture, RectangleSelfOr_Large)
    ->Arg(10000)->Arg(100000)->Iterations(3);

// Million-scale: single-shot stress tests
BENCHMARK_DEFINE_F(RectangleSelfOrFixture, RectangleSelfOr_Million)(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(self_or(poly));
}
BENCHMARK_REGISTER_F(RectangleSelfOrFixture, RectangleSelfOr_Million)
    ->Arg(500000)->Arg(1000000)->Iterations(1)->Unit(benchmark::kMillisecond);

// ============================================================================
// WKT union — two triangles
// ============================================================================
static void BM_Union_Triangles(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((-59 867,-36 492,-182 486,-59 867)))", a);
    bg::read_wkt("MULTIPOLYGON(((-220 877,-54 821,-402 541,-808 638,-220 877)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK(BM_Union_Triangles);

// ============================================================================
// WKT union — square with hole + overlapping square
// ============================================================================
static void BM_Union_Hole(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))", a);
    bg::read_wkt("MULTIPOLYGON(((2 2, 2 4, 4 4, 4 2, 2 2)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK(BM_Union_Hole);

// ============================================================================
// WKT union — two polygons with collinear edges
// ============================================================================
static void BM_Union_Collinear(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0, 1 1, 2 1, 2 2, 3 3, 3 0, 0 0)))", a);
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 3, 3 3, 2 2, 1 2, 1 1, 0 0)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK(BM_Union_Collinear);

// ============================================================================
// WKT intersection — two squares
// ============================================================================
static void BM_Intersection_Squares(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 2, 2 2, 2 0, 0 0)))", a);
    bg::read_wkt("MULTIPOLYGON(((1 1, 1 3, 3 3, 3 1, 1 1)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(intersection(a, b));
}
BENCHMARK(BM_Intersection_Squares);

// ============================================================================
// WKT intersection — disjoint (returns empty)
// ============================================================================
static void BM_Intersection_Disjoint(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))", a);
    bg::read_wkt("MULTIPOLYGON(((2 2, 2 3, 3 3, 3 2, 2 2)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(intersection(a, b));
}
BENCHMARK(BM_Intersection_Disjoint);

// ============================================================================
// WKT XOR — squares with hole
// ============================================================================
static void BM_Xor_Squares(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))", a);
    bg::read_wkt("MULTIPOLYGON(((2 2, 2 4, 4 4, 4 2, 2 2)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(xor_(a, b));
}
BENCHMARK(BM_Xor_Squares);

// ============================================================================
// WKT difference — large square minus middle
// ============================================================================
static void BM_Difference_Hole(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1), (6 6, 8 6, 8 8, 6 8, 6 6)))", a);
    bg::read_wkt("MULTIPOLYGON(((2 2, 2 7, 7 7, 7 2, 2 2)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(difference(a, b));
}
BENCHMARK(BM_Difference_Hole);

// ============================================================================
// Complex union — large WKT from test suite
// ============================================================================
static void BM_ComplexUnion(benchmark::State& state) {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((-1461 -786,-1417 -833,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1030 -1246,-1608 -939,-1461 -786),(-1058 -1133,-1244 -963,-1301 -1039,-1058 -1133)),((-578 -670,-688 -679,-802 -443,-826 -400,-578 -670)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))", a);
    bg::read_wkt("MULTIPOLYGON(((-1450 -1280, -1450 -800, -1200 -1000, -1000 -1280, -1450 -1280)))", b);
    for (auto _ : state)
        benchmark::DoNotOptimize(add(a, b));
}
BENCHMARK(BM_ComplexUnion);

BENCHMARK_MAIN();
