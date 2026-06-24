#include "best_clipper.hpp"
#include <benchmark/benchmark.h>

using namespace best_clipper;
#include <clipper2/clipper.h>

#include <algorithm>
#include <cmath>
#include <random>

// ============================================================================
// Random non-Manhattan 11-gons uniformly distributed in [0,1000]×[0,1000].
// Each polygon is a randomly perturbed regular 11-gon. Vertices are CW for
// Boost compatibility. Edges are slanted — no axis-aligned shortcuts.
// ============================================================================
namespace {

std::mt19937 rng(42);
constexpr int kSides = 66;
constexpr int32_t kExtent = 1000;

polygon gen_poly(int32_t cx, int32_t cy, int32_t radius, std::mt19937 &rng) {
  std::uniform_int_distribution<int32_t> jitter(-2, 2);
  polygon poly;
  auto &outer = poly.outer();
  for (int j = 0; j < kSides; j++) {
    double a = 2.0 * 3.141592653589793 * j / kSides;
    int32_t vx = cx + (int32_t)(radius * std::cos(a)) + jitter(rng);
    int32_t vy = cy + (int32_t)(radius * std::sin(a)) + jitter(rng);
    outer.push_back({vx, vy});
  }
  std::reverse(outer.begin(), outer.end()); // CW for Boost
  outer.push_back(outer.front());           // close
  return poly;
}

struct RandPolygons {
  multi_polygon polys;

  // Generate n random 11-gons uniformly in [0, kExtent]².
  // Radius 8..20, jitter ±2. Polygons may overlap at high n.
  explicit RandPolygons(int n) {
    polys.reserve(n);
    std::uniform_int_distribution<int32_t> center(60, kExtent - 60);
    std::uniform_int_distribution<int32_t> radius_dist(8, 20);
    for (int i = 0; i < n; i++)
      polys.push_back(
          gen_poly(center(rng), center(rng), radius_dist(rng), rng));
  }
};

// Non-overlapping variant: places polygons on a jittered grid so each
// multi_polygon is strictly valid (no within-set overlap).
struct RandPolygonsGrid {
  multi_polygon polys;

  // n polygons placed on a grid with spacing = 100, jitter ±15 in x/y.
  // Max polygon extent = 20+2 = 22, so min separation = 100-30-44 = 26 > 0.
  // Only works for n ≤ 100 (10×10 grid in 1000×1000).
  explicit RandPolygonsGrid(int n, uint32_t seed) {
    rng.seed(seed);
    polys.reserve(n);
    int cols = (int)std::ceil(std::sqrt((double)n));
    std::uniform_int_distribution<int32_t> jitter(-15, 15);
    std::uniform_int_distribution<int32_t> radius_dist(8, 20);
    for (int i = 0; i < n; i++) {
      int row = i / cols, col = i % cols;
      int32_t cx = 80 + col * 90 + jitter(rng);
      int32_t cy = 80 + row * 90 + jitter(rng);
      polys.push_back(gen_poly(cx, cy, radius_dist(rng), rng));
    }
  }
};

Clipper2Lib::Paths64 to_paths64(const multi_polygon &mp) {
  Clipper2Lib::Paths64 paths;
  for (const auto &poly : mp) {
    Clipper2Lib::Path64 outer;
    for (const auto &pt : poly.outer())
      outer.push_back({pt.x(), pt.y()});
    paths.push_back(std::move(outer));
    for (const auto &inner : poly.inners()) {
      Clipper2Lib::Path64 path;
      for (const auto &pt : inner)
        path.push_back({pt.x(), pt.y()});
      paths.push_back(std::move(path));
    }
  }
  return paths;
}

} // namespace

// ============================================================================
// Union benchmark
// ============================================================================
struct CmpUnionFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    RandPolygonsGrid ta(n, 42), tb(n, 99);
    a = std::move(ta.polys);
    b = std::move(tb.polys);
    paths_a = to_paths64(a);
    paths_b = to_paths64(b);
    cached_n = n;
  }
  multi_polygon a, b;
  Clipper2Lib::Paths64 paths_a, paths_b;
  int cached_n = 0;
};
#define CMP_UNION_ARGS ->Arg(10)->Arg(50)->Arg(100)->Arg(500)->Arg(1000)

BENCHMARK_DEFINE_F(CmpUnionFixture,
                   BestClipper_Union)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(union_(a, b));
}
BENCHMARK_REGISTER_F(CmpUnionFixture, BestClipper_Union) CMP_UNION_ARGS;

BENCHMARK_DEFINE_F(CmpUnionFixture, Clipper2_Union)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(
        Clipper2Lib::Union(paths_a, paths_b, Clipper2Lib::FillRule::NonZero));
}
BENCHMARK_REGISTER_F(CmpUnionFixture, Clipper2_Union) CMP_UNION_ARGS;

// ============================================================================
// Intersection benchmark
// ============================================================================
struct CmpIntersectionFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    RandPolygonsGrid ta(n, 42), tb(n, 99);
    a = std::move(ta.polys);
    b = std::move(tb.polys);
    paths_a = to_paths64(a);
    paths_b = to_paths64(b);
    cached_n = n;
  }
  multi_polygon a, b;
  Clipper2Lib::Paths64 paths_a, paths_b;
  int cached_n = 0;
};
#define CMP_INTER_ARGS ->Arg(10)->Arg(50)->Arg(100)->Arg(500)->Arg(1000)
BENCHMARK_DEFINE_F(CmpIntersectionFixture,
                   BestClipper_Intersection)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(intersection(a, b));
}
BENCHMARK_REGISTER_F(CmpIntersectionFixture, BestClipper_Intersection)
CMP_INTER_ARGS;
BENCHMARK_DEFINE_F(CmpIntersectionFixture,
                   Clipper2_Intersection)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(Clipper2Lib::Intersect(
        paths_a, paths_b, Clipper2Lib::FillRule::NonZero));
}
BENCHMARK_REGISTER_F(CmpIntersectionFixture, Clipper2_Intersection)
CMP_INTER_ARGS;

// ============================================================================
// XOR benchmark
// ============================================================================
struct CmpXorFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    RandPolygonsGrid ta(n, 42), tb(n, 99);
    a = std::move(ta.polys);
    b = std::move(tb.polys);
    paths_a = to_paths64(a);
    paths_b = to_paths64(b);
    cached_n = n;
  }
  multi_polygon a, b;
  Clipper2Lib::Paths64 paths_a, paths_b;
  int cached_n = 0;
};
#define CMP_XOR_ARGS ->Arg(10)->Arg(50)->Arg(100)->Arg(500)->Arg(1000)
BENCHMARK_DEFINE_F(CmpXorFixture, BestClipper_Xor)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(symmetric_difference(a, b));
}
BENCHMARK_REGISTER_F(CmpXorFixture, BestClipper_Xor) CMP_XOR_ARGS;
BENCHMARK_DEFINE_F(CmpXorFixture, Clipper2_Xor)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(
        Clipper2Lib::Xor(paths_a, paths_b, Clipper2Lib::FillRule::NonZero));
}
BENCHMARK_REGISTER_F(CmpXorFixture, Clipper2_Xor) CMP_XOR_ARGS;

// ============================================================================
// Difference benchmark
// ============================================================================
struct CmpDifferenceFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    RandPolygonsGrid ta(n, 42), tb(n, 99);
    a = std::move(ta.polys);
    b = std::move(tb.polys);
    paths_a = to_paths64(a);
    paths_b = to_paths64(b);
    cached_n = n;
  }
  multi_polygon a, b;
  Clipper2Lib::Paths64 paths_a, paths_b;
  int cached_n = 0;
};
#define CMP_DIFF_ARGS ->Arg(10)->Arg(50)->Arg(100)->Arg(500)->Arg(1000)
BENCHMARK_DEFINE_F(CmpDifferenceFixture,
                   BestClipper_Difference)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(difference(a, b));
}
BENCHMARK_REGISTER_F(CmpDifferenceFixture, BestClipper_Difference)
CMP_DIFF_ARGS;
BENCHMARK_DEFINE_F(CmpDifferenceFixture,
                   Clipper2_Difference)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(Clipper2Lib::Difference(
        paths_a, paths_b, Clipper2Lib::FillRule::NonZero));
}
BENCHMARK_REGISTER_F(CmpDifferenceFixture, Clipper2_Difference) CMP_DIFF_ARGS;

// ============================================================================
// Self-union benchmark — single set
// ============================================================================
struct CmpSelfOrFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    RandPolygonsGrid t(n, 42);
    poly = std::move(t.polys);
    paths = to_paths64(poly);
    cached_n = n;
  }
  multi_polygon poly;
  Clipper2Lib::Paths64 paths;
  int cached_n = 0;
};
#define CMP_SELFOR_ARGS ->Arg(10)->Arg(50)->Arg(100)->Arg(500)->Arg(1000)
BENCHMARK_DEFINE_F(CmpSelfOrFixture,
                   BestClipper_SelfOr)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(robust_self_or(poly));
}
BENCHMARK_REGISTER_F(CmpSelfOrFixture, BestClipper_SelfOr) CMP_SELFOR_ARGS;
BENCHMARK_DEFINE_F(CmpSelfOrFixture, Clipper2_SelfOr)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(
        Clipper2Lib::Union(paths, Clipper2Lib::FillRule::Positive));
}
BENCHMARK_REGISTER_F(CmpSelfOrFixture, Clipper2_SelfOr) CMP_SELFOR_ARGS;

BENCHMARK_MAIN();
