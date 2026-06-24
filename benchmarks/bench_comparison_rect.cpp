#include "best_clipper.hpp"
#include <benchmark/benchmark.h>
#include <clipper2/clipper.h>
#include <random>

using namespace best_clipper;

// ============================================================================
// Uniform random rectangles in [0, 1000]^2 — realistic spatial distribution
// ============================================================================
static polygon make_rect(int32_t x, int32_t y, int32_t w, int32_t h) {
  polygon p;
  p.outer() = {{x, y}, {x, y + h}, {x + w, y + h}, {x + w, y}, {x, y}};
  return p;
}

static Clipper2Lib::Paths64 to_paths64(const multi_polygon &mp) {
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

namespace {
std::mt19937 rng(42);
}

// ============================================================================
// Union: a and b are independent random rect sets, uniformly in [0,1000]^2
// ============================================================================
struct RectUnionFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    a.clear();
    b.clear();
    a.reserve(n);
    b.reserve(n);
    std::uniform_int_distribution<int32_t> pos(0, (int32_t)(7.0 * std::sqrt(n)));
    std::uniform_int_distribution<int32_t> sz(2, 4);
    for (int i = 0; i < n; i++)
      a.push_back(make_rect(pos(rng), pos(rng), sz(rng), sz(rng)));
    for (int i = 0; i < n; i++)
      b.push_back(make_rect(pos(rng), pos(rng), sz(rng), sz(rng)));
    paths_a = to_paths64(a);
    paths_b = to_paths64(b);
    cached_n = n;
  }
  multi_polygon a, b;
  Clipper2Lib::Paths64 paths_a, paths_b;
  int cached_n = 0;
};

BENCHMARK_DEFINE_F(RectUnionFixture, BestClipper)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(union_(a, b));
}
BENCHMARK_REGISTER_F(RectUnionFixture, BestClipper)
    ->Arg(10)
    ->Arg(100)
    ->Arg(1000)
    ->Arg(10000)
    ->Arg(100000);

BENCHMARK_DEFINE_F(RectUnionFixture, Clipper2)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(
        Clipper2Lib::Union(paths_a, paths_b, Clipper2Lib::FillRule::NonZero));
}
BENCHMARK_REGISTER_F(RectUnionFixture, Clipper2)
    ->Arg(10)
    ->Arg(100)
    ->Arg(1000)
    ->Arg(10000)
    ->Arg(100000);

// ============================================================================
// Self-union: N random rectangles, uniform in plane
// ============================================================================
struct RectSelfOrFixture : benchmark::Fixture {
  void SetUp(const benchmark::State &state) override {
    int n = state.range(0);
    if (cached_n == n)
      return;
    poly.clear();
    poly.reserve(n);
    std::uniform_int_distribution<int32_t> pos(0, (int32_t)(7.0 * std::sqrt(n)));
    std::uniform_int_distribution<int32_t> sz(2, 4);
    for (int i = 0; i < n; i++)
      poly.push_back(make_rect(pos(rng), pos(rng), sz(rng), sz(rng)));
    paths = to_paths64(poly);
    cached_n = n;
  }
  multi_polygon poly;
  Clipper2Lib::Paths64 paths;
  int cached_n = 0;
};

BENCHMARK_DEFINE_F(RectSelfOrFixture, BestClipper)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(robust_self_or(poly));
}
BENCHMARK_REGISTER_F(RectSelfOrFixture, BestClipper)
    ->Arg(10)
    ->Arg(100)
    ->Arg(1000)
    ->Arg(10000)
    ->Arg(100000)
    ->Arg(500000);

BENCHMARK_DEFINE_F(RectSelfOrFixture, Clipper2)(benchmark::State &state) {
  for (auto _ : state)
    benchmark::DoNotOptimize(
        Clipper2Lib::Union(paths, Clipper2Lib::FillRule::Positive));
}
BENCHMARK_REGISTER_F(RectSelfOrFixture, Clipper2)
    ->Arg(10)
    ->Arg(100)
    ->Arg(1000)
    ->Arg(10000)
    ->Arg(100000)
    ->Arg(500000);

BENCHMARK_MAIN();
