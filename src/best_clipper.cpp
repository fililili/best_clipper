#include "include/best_clipper.hpp"
#include "filter_half_chains.hpp"
#include "include/chain_builder.hpp"
#include "include/connect_half_chains.hpp"
#include "include/half_chain_graph.hpp"
#include "include/output_builder.hpp"

#include <chrono>
#include <cstdio>
#include <utility>

namespace best_clipper {

// ---------------------------------------------------------------------------
// Append helpers - each ring writes points then pushes points.size().
// Ring i spans points[offsets[i] .. offsets[i+1]), offsets starts with {0}.
// ---------------------------------------------------------------------------

static void append_ring(std::vector<point> &points,
                        std::vector<std::size_t> &offsets, const ring &r) {
  if (r.empty())
    return;
  for (const auto &pt : r)
    points.push_back(pt);
  if (!bg::equals(r[0], r.back()))
    points.push_back(r[0]);
  offsets.push_back(points.size());
}

static void append_ring_rev(std::vector<point> &points,
                            std::vector<std::size_t> &offsets, const ring &r) {
  if (r.empty())
    return;
  bool ring_closed = bg::equals(r[0], r.back());
  points.push_back(r[0]);
  int last_original = ring_closed ? (int)r.size() - 2 : (int)r.size() - 1;
  for (int i = last_original; i > 0; i--)
    points.push_back(r[i]);
  points.push_back(r[0]);
  offsets.push_back(points.size());
}

static void append_polygon(std::vector<point> &points,
                           std::vector<std::size_t> &offsets,
                           const polygon &poly) {
  append_ring(points, offsets, poly.outer());
  for (const auto &inner : poly.inners())
    append_ring(points, offsets, inner);
}

static void append_polygon_rev(std::vector<point> &points,
                               std::vector<std::size_t> &offsets,
                               const polygon &poly) {
  append_ring_rev(points, offsets, poly.outer());
  for (const auto &inner : poly.inners())
    append_ring_rev(points, offsets, inner);
}

static void append_multi_polygon(std::vector<point> &points,
                                 std::vector<std::size_t> &offsets,
                                 const multi_polygon &mp) {
  for (const auto &poly : mp)
    append_polygon(points, offsets, poly);
}

static void append_multi_polygon_rev(std::vector<point> &points,
                                     std::vector<std::size_t> &offsets,
                                     const multi_polygon &mp) {
  for (const auto &poly : mp)
    append_polygon_rev(points, offsets, poly);
}

// ---------------------------------------------------------------------------
// collect_segments - gather rings into points + [begin, end) offsets
// ---------------------------------------------------------------------------

static auto collect_segments(const multi_polygon &mp) {
  std::vector<point> points;
  std::vector<std::size_t> offsets = {0};
  append_multi_polygon(points, offsets, mp);
  return std::pair{std::move(points), std::move(offsets)};
}

static auto collect_segments(const multi_polygon &mp1,
                             const multi_polygon &mp2) {
  auto [points, offsets] = collect_segments(mp1);
  append_multi_polygon(points, offsets, mp2);
  return std::pair{std::move(points), std::move(offsets)};
}

inline auto run_pipeline(std::vector<point> points,
                         std::vector<std::size_t> offsets, auto filter) {
  auto [chains, out_offsets, out_chains, in_offsets, in_chains] =
      build_chains_from_input(points, offsets);

  auto [bucket_half_chains, bucket_half_chains_offsets] =
      build_half_chain_graph(chains, out_offsets, out_chains, in_offsets,
                             in_chains);

  auto half_chain_relations =
      build_half_chain_relations(chains, std::move(bucket_half_chains),
                                 std::move(bucket_half_chains_offsets));

  auto winding = compute_winding(chains, half_chain_relations);

  auto survive = filter_survive(winding, filter);

  auto result = build_output(chains, half_chain_relations, survive);
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
  append_multi_polygon_rev(points, offsets, b);
  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w > 0; });
}

multi_polygon robust_self_or(const multi_polygon &a) {
  auto [points, offsets] = collect_segments(a);
  return run_pipeline(std::move(points), std::move(offsets),
                      [](int w) { return w > 0; });
}

} // namespace best_clipper
