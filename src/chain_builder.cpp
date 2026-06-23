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
  std::vector<std::size_t> seg_idx;
  seg_boxes.reserve(points.size() - 1);
  seg_idx.reserve(points.size() - 1);
  {
    for (std::size_t ri = 0; ri + 1 < offsets.size(); ri++) {
      std::size_t rb = offsets[ri];
      std::size_t re = offsets[ri + 1];
      for (std::size_t i = rb; i + 1 < re; i++) {
        seg_boxes.push_back(
            bg::return_envelope<box>(segment{points[i], points[i + 1]}));
        seg_idx.push_back(i);
      }
      if (ri + 2 < offsets.size()) {
        seg_boxes.push_back(box{point{1, 1}, point{0, 0}});
        seg_idx.push_back(rb);
      }
    }
  }
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
        segments_box_grid.query_intersects(box_i, [&](std::size_t pos) {
          std::size_t j = seg_idx[pos];
          if (i <= j + 1)
            return; // adjacent edges share a vertex, cannot produce a new
                    // intersection
          if (auto p = get_intersection(si, segment{points[j], points[j + 1]}))
            hot_pixels.push_back(p.value());
          // todo: if no more intersections are possible, we could return chains from input direclty
        });
      }
    }
  }

  // Dedup hot_pixels: group by grid cell, sort within cell, then unique.
  {
    auto &g = segments_box_grid;
    std::size_t num_cells = (std::size_t)(g._x_cells * g._y_cells);
    std::vector<std::size_t> pixel_cell(hot_pixels.size());
    for (std::size_t pi = 0; pi < hot_pixels.size(); pi++) {
      auto &p = hot_pixels[pi];
      pixel_cell[pi] =
          (std::size_t)(((bg::get<1>(p) - g._min_y) / g._cell_size) *
                            g._x_cells +
                        (bg::get<0>(p) - g._min_x) / g._cell_size);
    }
    std::vector<std::size_t> cell_counts(num_cells, 0);
    for (auto c : pixel_cell)
      cell_counts[c]++;
    std::vector<std::size_t> begins(num_cells + 1, 0);
    for (std::size_t c = 0; c < num_cells; c++)
      begins[c + 1] = begins[c] + cell_counts[c];
    std::vector<std::size_t> cell_items(hot_pixels.size());
    auto cursors = begins;
    for (std::size_t pi = 0; pi < hot_pixels.size(); pi++)
      cell_items[cursors[pixel_cell[pi]]++] = pi;
    std::vector<point> unique;
    unique.reserve(hot_pixels.size());
    for (std::size_t c = 0; c < num_cells; c++) {
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
    auto min_corner = point{x > INT32_MIN ? x - 1 : INT32_MIN,
                            y > INT32_MIN ? y - 1 : INT32_MIN};
    auto max_corner = point{x < INT32_MAX ? x + 1 : INT32_MAX,
                            y < INT32_MAX ? y + 1 : INT32_MAX};
    candidates.clear();
    segments_box_grid.query_intersects(box{min_corner, max_corner},
                                       [&](std::size_t pos) {
                                         std::size_t seg_start = seg_idx[pos];
                                         if (last_seen[seg_start] == pi)
                                           return;
                                         last_seen[seg_start] = pi;
                                         candidates.push_back(seg_start);
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
        std::sort(std::begin(pixels) + segs[i], std::begin(pixels) + segs[i + 1], [&](auto pi, auto pj) {
          return less_by_segment{seg}(hot_pixels[pi], hot_pixels[pj]);
        });
        for(auto k = segs[i]; k + 1 < segs[i + 1]; ++k) {
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

std::tuple<std::vector<point>, chain_group>
build_chains_from_input(const std::vector<point> &points,
                        const std::vector<std::size_t> &offsets) {
#ifndef NDEBUG
  assert(offsets.size() >= 2 && offsets[0] == 0 &&
         offsets.back() == points.size());
  for (std::size_t ri = 0; ri + 1 < offsets.size(); ri++)
    assert(bg::equals(points[offsets[ri]], points[offsets[ri + 1] - 1]));
#endif
  auto t0 = std::chrono::high_resolution_clock::now();
  auto [edges, hot_pixels] = construct_graph(points, offsets);
  auto t1 = std::chrono::high_resolution_clock::now();

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
  auto t2 = std::chrono::high_resolution_clock::now();

  auto chains = build_chains(sorted_edges, hot_pixels.size());

  auto t3 = std::chrono::high_resolution_clock::now();
  auto ms = [](auto d) {
    return std::chrono::duration<double, std::milli>(d).count();
  };
  std::fprintf(stderr,
               "  [edges] graph=%.1fms dedup=%.1fms chains=%.1fms (hp=%zu)\n",
               ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), hot_pixels.size());

  return std::tuple{std::move(hot_pixels), std::move(chains)};
}

} // namespace best_clipper
