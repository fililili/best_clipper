#include "include/chain_builder.hpp"
#include "include/graph_helper.hpp"
#include "include/snap_rounding_helper.hpp"
#include "include/uniform_grid.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

namespace best_clipper {

namespace bg = boost::geometry;

// ---------------------------------------------------------------------------
// construct_graph — returns edges and hot_pixels
// ---------------------------------------------------------------------------

std::tuple<std::vector<edge_t>, std::vector<point>>
construct_graph(const std::vector<point> &points,
                const std::vector<std::size_t> &offsets) {
  // Build grid: point index i is the segment (points[i], points[i+1]).
  // offsets[i] = start of ring i; ring i spans [offsets[i], offsets[i+1]).
  // offsets.back() = points.size() is a sentinel.  Barrier segments between
  // rings get an invalid bbox so the grid skips them naturally.
  std::vector<box> seg_boxes;
  seg_boxes.reserve(points.size());
  {
    for (std::size_t ri = 0; ri + 1 < offsets.size(); ri++) {
      std::size_t rb = offsets[ri];
      std::size_t re = offsets[ri + 1];
      for (std::size_t i = rb; i + 1 < re; i++) {
        seg_boxes.push_back(
            bg::return_envelope<box>(segment{points[i], points[i + 1]}));
      }
      seg_boxes.push_back(box{point{1, 1}, point{0, 0}});
    }
  }
  assert(seg_boxes.size() == points.size());
  best_clipper::uniform_grid::grid segments_box_grid(seg_boxes);

  // Hot pixels from all points
  std::vector<point> hot_pixels = points;

  // Find intersections between segments (skip offset boundaries)
  {
    for (std::size_t ri = 0; ri + 1 < offsets.size(); ri++) {
      std::size_t rb = offsets[ri];
      std::size_t re = offsets[ri + 1];
      for (std::size_t i = rb; i + 1 < re; i++) {
        auto si = segment{points[i], points[i + 1]};
        box box_i = bg::return_envelope<box>(si);
        segments_box_grid.query_intersects(box_i, [&](std::size_t j) {
          if (i <= j + 1) {
            return; // adjacent edges share a vertex, cannot produce a new
                    // intersection
          }
          if (!bbox_overlap(box_i, seg_boxes[j])) {
            return; // grid cell false positive
          }
          if (auto p = get_intersection(si, segment{points[j], points[j + 1]}))
            hot_pixels.push_back(p.value());
          // todo: if no more intersections are possible, we could return chains
          // from input direclty
        });
      }
    }
  }

  // Dedup hot_pixels: group by grid cell, sort within cell, then unique.
  {
    std::vector<std::size_t> pixel_cell(hot_pixels.size());
    for (std::size_t pi = 0; pi < hot_pixels.size(); pi++) {
      pixel_cell[pi] = segments_box_grid.get_flat_index(
          segments_box_grid.cell_x(bg::get<0>(hot_pixels[pi])),
          segments_box_grid.cell_y(bg::get<1>(hot_pixels[pi])));
    }
    std::vector<std::size_t> cell_counts(segments_box_grid.cells(), 0);
    for (auto c : pixel_cell)
      cell_counts[c]++;
    std::vector<std::size_t> begins(segments_box_grid.cells() + 1, 0);
    for (std::size_t c = 0; c < segments_box_grid.cells(); c++)
      begins[c + 1] = begins[c] + cell_counts[c];
    std::vector<std::size_t> cell_items(hot_pixels.size());
    auto cursors = begins;
    for (std::size_t pi = 0; pi < hot_pixels.size(); pi++)
      cell_items[cursors[pixel_cell[pi]]++] = pi;
    std::vector<point> unique;
    unique.reserve(hot_pixels.size());
    for (std::size_t c = 0; c < segments_box_grid.cells(); c++) {
      auto b = begins[c], e = begins[c + 1];
      std::sort(cell_items.begin() + b, cell_items.begin() + e,
                [&](std::size_t a, std::size_t b) {
                  return std::pair{bg::get<0>(hot_pixels[a]),
                                   bg::get<1>(hot_pixels[a])} <
                         std::pair{bg::get<0>(hot_pixels[b]),
                                   bg::get<1>(hot_pixels[b])};
                });
      for (auto it = cell_items.begin() + b; it != cell_items.begin() + e;
           ++it) {
        if (it == cell_items.begin() + b ||
            !bg::equals(hot_pixels[*it], hot_pixels[*(it - 1)]))
          unique.push_back(hot_pixels[*it]);
      }
    }
    hot_pixels = std::move(unique);
  }

  // Build segment-pixel pairs.  A hot pixel is always on ≥2 segments
  // (original vertex or intersection), so when the grid returns exactly 2
  // candidates we skip the expensive is_point_on_segment check.
  std::vector<std::pair<std::size_t, std::size_t>> segment_pixel_pairs;
  std::vector<std::size_t> last_seen(points.size(), ~0ULL);
  std::vector<std::size_t> candidates;
  for (std::size_t pi = 0; pi < hot_pixels.size(); pi++) {
    int32_t x = bg::get<0>(hot_pixels[pi]), y = bg::get<1>(hot_pixels[pi]);
    auto query_box = box{{x - 1, y - 1}, {x + 1, y + 1}};
    candidates.clear();
    segments_box_grid.query_intersects(query_box, [&](std::size_t pos) {
      if (last_seen[pos] == pi)
        return;
      last_seen[pos] = pi;
      if (!bbox_overlap(query_box, seg_boxes[pos]))
        return;
      candidates.push_back(pos);
    });
    if (candidates.size() == 2) {
      segment_pixel_pairs.emplace_back(candidates[0], pi);
      segment_pixel_pairs.emplace_back(candidates[1], pi);
    } else {
      for (auto seg_start : candidates) {
        if (is_point_on_segment(hot_pixels[pi], segment{points[seg_start],
                                                        points[seg_start + 1]}))
          segment_pixel_pairs.emplace_back(seg_start, pi);
      }
    }
  }

  // Build edges: sort pixels along each segment, connect adjacent pairs
  std::vector<edge_t> edges;
  {
    auto [segs, pixels] = bucket_sort(
        segment_pixel_pairs, points.size(), [](auto val) { return val.first; },
        [](auto val) { return val.second; });
    for (std::size_t ri = 0; ri + 1 < offsets.size(); ri++) {
      std::size_t rb = offsets[ri];
      std::size_t re = offsets[ri + 1];
      for (std::size_t i = rb; i + 1 < re; i++) {
        auto seg = segment{points[i], points[i + 1]};
        std::sort(std::begin(pixels) + segs[i],
                  std::begin(pixels) + segs[i + 1], [&](auto pi, auto pj) {
                    return less_by_segment{seg}(hot_pixels[pi], hot_pixels[pj]);
                  });
        for (auto k = segs[i]; k + 1 < segs[i + 1]; ++k) {
          edges.emplace_back(pixels[k], pixels[k + 1]);
        }
      }
    }
  }
  return {std::move(edges), std::move(hot_pixels)};
}

// ---------------------------------------------------------------------------
// edges_to_power
// ---------------------------------------------------------------------------

std::vector<edge_with_power_t> edges_to_power(std::vector<edge_t> edges) {
  std::vector<edge_with_power_t> r(edges.size() * 2);
  for (std::size_t i = 0; i < edges.size(); i++) {
    auto s = edges[i].start, e = edges[i].end;
    r[2 * i] = {s, e, 1};
    r[2 * i + 1] = {e, s, -1};
  }
  return r;
}

// ---------------------------------------------------------------------------
// unique_edges
// ---------------------------------------------------------------------------

std::vector<edge_with_power_t>
unique_edges(std::vector<edge_with_power_t> edges, std::size_t num_vertices) {
  auto [locs, ordered] = bucket_sort(
      std::move(edges), num_vertices,
      [](const edge_with_power_t &e) { return e.start; },
      [](const edge_with_power_t &e) { return e; });
  std::vector<edge_with_power_t> result;
  for (std::size_t i = 0; i < num_vertices; i++) {
    auto current_begin = std::begin(ordered) + locs[i],
         current_end = std::begin(ordered) + locs[i + 1];
    if (current_begin == current_end)
      continue;
    std::sort(current_begin, current_end,
              [](const edge_with_power_t &a, const edge_with_power_t &b) {
                return a.end < b.end;
              });
    for (auto cur = current_begin; cur != current_end;) {
      int sum = 0;
      auto next = cur;
      while (next != current_end && next->end == cur->end) {
        sum += next->power;
        ++next;
      }
      if (sum > 0)
        result.push_back({cur->start, cur->end, sum});
      cur = next;
    }
  }
  return result;
}

// ---------------------------------------------------------------------------
// build_chains_from_input
// ---------------------------------------------------------------------------

std::tuple<chain_group, std::vector<std::size_t>, std::vector<std::size_t>,
           std::vector<std::size_t>, std::vector<std::size_t>>
build_chains_from_input(const std::vector<point> &points,
                        const std::vector<std::size_t> &offsets) {
  assert(offsets.size() >= 1);
  assert(offsets[0] == 0);
  assert(offsets.back() == points.size());
#ifndef NDEBUG
  for (std::size_t ri = 0; ri + 1 < offsets.size(); ri++)
    assert(bg::equals(points[offsets[ri]], points[offsets[ri + 1] - 1]));
#endif
  auto [edges, hot_pixels] = construct_graph(points, offsets);

  auto sorted_edges =
      unique_edges(edges_to_power(std::move(edges)), hot_pixels.size());

#ifndef NDEBUG
  {
    std::vector<int> balance(hot_pixels.size());
    for (const auto &e : sorted_edges)
      balance[e.start] += e.power, balance[e.end] -= e.power;
    for (std::size_t v = 0; v < hot_pixels.size(); v++)
      assert(balance[v] == 0);
  }
#endif
  auto [chain_indexes, chain_offsets, chain_powers, chain_out_offsets,
        chain_out_chains, chain_in_offsets, chain_in_chains] =
      build_chains(sorted_edges, hot_pixels.size());
  return {chain_group{std::move(hot_pixels), std::move(chain_indexes),
                      std::move(chain_offsets), std::move(chain_powers)},
          std::move(chain_out_offsets), std::move(chain_out_chains),
          std::move(chain_in_offsets), std::move(chain_in_chains)};
}

} // namespace best_clipper
