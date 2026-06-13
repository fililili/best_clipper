#include "include/best_clipper.hpp"
#include "include/chain_builder.hpp"
#include "include/exterior_winding.hpp"
#include "include/half_chain_graph.hpp"
#include "include/output_builder.hpp"

#include <chrono>
#include <cstdio>
#include <utility>

namespace best_clipper {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static auto collect_segments(const multi_polygon &mp) {
  std::vector<point> points;
  std::vector<std::size_t> offsets;

  auto add_ring = [&](const auto &ring) {
    if (ring.empty())
      return;
    if (!points.empty())
      offsets.push_back(points.size() - 1);
    for (const auto &pt : ring)
      points.push_back(pt);
  };

  for (const auto &poly : mp) {
    add_ring(poly.outer());
    for (const auto &inner : poly.inners())
      add_ring(inner);
  }

  return std::pair{std::move(points), std::move(offsets)};
}

static auto collect_segments(const multi_polygon &mp1,
                             const multi_polygon &mp2) {
  auto [points, offsets] = collect_segments(mp1);
  auto add_mp = [&](const auto &mp) {
    for (const auto &poly : mp) {
      auto add_ring = [&](const auto &ring) {
        if (ring.empty())
          return;
        if (!points.empty())
          offsets.push_back(points.size() - 1);
        for (const auto &pt : ring)
          points.push_back(pt);
      };
      add_ring(poly.outer());
      for (const auto &inner : poly.inners())
        add_ring(inner);
    }
  };
  add_mp(mp2);
  return std::pair{std::move(points), std::move(offsets)};
}

inline std::vector<bool> filter_survive(const std::vector<int> &winding,
                                        auto filter_fn) {
  std::vector<bool> survive(winding.size());
  for (std::size_t i = 0; i < winding.size(); i++)
    survive[i] = filter_fn(winding[i]);
  return survive;
}

inline auto run_pipeline(std::vector<point> points,
                         std::vector<std::size_t> offsets, auto filter) {
  using clock = std::chrono::high_resolution_clock;
  auto t0 = clock::now();

  auto [hot_pixels, chains] = build_chains_from_input(points, offsets);
  auto t1 = clock::now();

  auto [sorted_half_chains, half_chain_begin, half_chain_end, next_half_chain,
        coplanar] = build_half_chain_graph(chains, hot_pixels);
  auto t2 = clock::now();

  auto [exterior_half_chains, ray_pairs] = find_exterior(
      chains, hot_pixels, sorted_half_chains, half_chain_begin, half_chain_end);
  auto t3 = clock::now();

  auto winding =
      compute_winding(chains, coplanar, ray_pairs, exterior_half_chains);
  auto t4 = clock::now();

  auto survive = filter_survive(winding, filter);

  auto result = build_output(chains, hot_pixels, std::move(next_half_chain),
                             std::move(survive), coplanar, ray_pairs);
  auto t5 = clock::now();

  auto ms = [](auto d) {
    return std::chrono::duration<double, std::milli>(d).count();
  };
  std::fprintf(stderr,
               "[pipeline] edges=%.1fms half_graph=%.1fms exterior=%.1fms "
               "winding=%.1fms output=%.1fms total=%.1fms\n",
               ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), ms(t4 - t3), ms(t5 - t4),
               ms(t5 - t0));
  return result;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

multi_polygon union_(const multi_polygon &a, const multi_polygon &b) {
  auto [points, offsets] = collect_segments(a, b);
  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w > 0; });
}

multi_polygon intersection(const multi_polygon &a, const multi_polygon &b) {
  auto [points, offsets] = collect_segments(a, b);
  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w > 1; });
}

multi_polygon symmetric_difference(const multi_polygon &a,
                                   const multi_polygon &b) {
  auto [points, offsets] = collect_segments(a, b);
  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w == 1; });
}

multi_polygon difference(const multi_polygon &a, const multi_polygon &b) {
  auto [points, offsets] = collect_segments(a);

  // Add b's rings in reverse order
  auto add_mp_rev = [&](const auto &mp) {
    for (const auto &poly : mp) {
      auto add_ring_rev = [&](const auto &ring) {
        if (ring.empty())
          return;
        if (!points.empty())
          offsets.push_back(points.size() - 1);
        points.push_back(ring[0]);
        for (int i = (int)ring.size() - 2; i >= 0; i--)
          points.push_back(ring[i]);
      };
      add_ring_rev(poly.outer());
      for (const auto &inner : poly.inners())
        add_ring_rev(inner);
    }
  };
  add_mp_rev(b);

  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w > 0; });
}

multi_polygon self_or(const multi_polygon &a) {
  auto [points, offsets] = collect_segments(a);
  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w > 0; });
}

} // namespace best_clipper
