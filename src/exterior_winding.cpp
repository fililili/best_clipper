#include "include/exterior_winding.hpp"
#include "include/graph_helper.hpp"
#include "include/snap_rounding_helper.hpp"
#include "include/uniform_grid.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

namespace best_clipper {

namespace bg = boost::geometry;

// ---------------------------------------------------------------------------
// chain_seg (internal to ray casting)
// ---------------------------------------------------------------------------

struct chain_seg {
  int32_t x1, y1, x2, y2;
  std::size_t half_chain_id;
};

// ---------------------------------------------------------------------------
// cast_ray_minus_x
// ---------------------------------------------------------------------------

std::size_t cast_ray_minus_x(coordinate_type vx, coordinate_type ray_y,
                             const best_clipper::uniform_grid::grid &seg_grid,
                             const std::vector<chain_seg> &seg_data) {

  bool any_hit = false;
  std::size_t hit_half_chain_id = (std::size_t)-1;
  coordinate_type best_hit_x1 = 0;
  coordinate_type best_hit_x2 = 0;
  coordinate_type best_hit_y1 = 0;
  coordinate_type best_hit_y2 = 0;
  auto point_on_seg_right = [](coordinate_type x, coordinate_type y,
                               coordinate_type x1, coordinate_type y1,
                               coordinate_type x2, coordinate_type y2) {
    return (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) < 0;
  };
  auto none_intersect_segs_compare_by_xray =
      [point_on_seg_right](coordinate_type x11, coordinate_type x12,
                           coordinate_type y11, coordinate_type y12,
                           coordinate_type x21, coordinate_type x22,
                           coordinate_type y21, coordinate_type y22) {
        if (y11 <= y21) {
          return point_on_seg_right(x21, y21, x11, y11, x12, y12);
        } else {
          return !point_on_seg_right(x11, y11, x21, y21, x22, y22);
        }
      };

  auto query_box =
      box{point{std::numeric_limits<coordinate_type>::min(), ray_y},
          point{vx, ray_y}};
  seg_grid.query_intersects(query_box, [&](std::size_t idx) {
    auto &seg = seg_data[idx];
    assert(seg.y1 < seg.y2);
    if (seg.y1 <= ray_y && ray_y < seg.y2 &&
        point_on_seg_right(vx, ray_y, seg.x1, seg.y1, seg.x2, seg.y2)) {
      if (!any_hit) {
        best_hit_x1 = seg.x1;
        best_hit_y1 = seg.y1;
        best_hit_x2 = seg.x2;
        best_hit_y2 = seg.y2;
        hit_half_chain_id = seg.half_chain_id;
        any_hit = true;
      } else if (none_intersect_segs_compare_by_xray(
                     best_hit_x1, best_hit_x2, best_hit_y1, best_hit_y2, seg.x1,
                     seg.x2, seg.y1, seg.y2)) {
        best_hit_x1 = seg.x1;
        best_hit_y1 = seg.y1;
        best_hit_x2 = seg.x2;
        best_hit_y2 = seg.y2;
        hit_half_chain_id = seg.half_chain_id;
      }
    }
  });

  return hit_half_chain_id;
}

// ---------------------------------------------------------------------------
// find_exterior
// ---------------------------------------------------------------------------

fe_tuple
find_exterior(const chain_build_result &chains,
              const std::vector<point> &hot_pixels,
              const std::vector<half_chain> &sorted_half_chains,
              const std::vector<std::size_t> &sorted_half_chains_offsets) {
  std::size_t num_vertices = hot_pixels.size();

  std::vector<std::pair<std::size_t, std::size_t>> vertex_edges;
  for (std::size_t c = 0; c + 1 < chains.offsets.size(); c++) {
    auto chain_begin_idx = chains.offsets[c],
         chain_end_idx = chains.offsets[c + 1];
    for (std::size_t k = chain_begin_idx; k + 1 < chain_end_idx; k++)
      vertex_edges.emplace_back(chains.indices[k], chains.indices[k + 1]);
  }

  auto component_ids = connected_components(num_vertices, vertex_edges);

  std::size_t num_components = 0;
  for (auto c : component_ids)
    num_components = std::max(num_components, c + 1);

  std::vector<size_t> leftmost_vertexes(num_components);
  std::vector<coordinate_type> min_xs(num_components,
                                      std::numeric_limits<int32_t>::max());
  for (std::size_t vertex = 0; vertex < component_ids.size(); ++vertex) {
    auto x = bg::get<0>(hot_pixels[vertex]);
    auto component_id = component_ids[vertex];

    if (!(sorted_half_chains_offsets[vertex] <
          sorted_half_chains_offsets[vertex + 1])) {
      continue;
    }
    if (x < min_xs[component_id]) {
      min_xs[component_id] = x;
      leftmost_vertexes[component_id] = vertex;
    }
  }

  std::vector<std::pair<std::size_t, std::size_t>> chain_edges;
  for (std::size_t c = 0; c + 1 < chains.offsets.size(); c++) {
    auto chain_begin_idx = chains.offsets[c],
         chain_end_idx = chains.offsets[c + 1];
    for (std::size_t k = chain_begin_idx; k + 1 < chain_end_idx; k++)
      vertex_edges.emplace_back(chains.indices[k], chains.indices[k + 1]);
  }

  std::vector<chain_seg> seg_data;
  std::vector<box> seg_boxes;
  size_t total_edges = 0;
  for (size_t ci = 0; ci + 1 < chains.offsets.size(); ci++)
    total_edges += chains.offsets[ci + 1] - chains.offsets[ci] - 1;
  seg_data.reserve(total_edges);
  seg_boxes.reserve(total_edges);
  for (size_t ci = 0; ci + 1 < chains.offsets.size(); ci++) {
    auto cbi = chains.offsets[ci], cei = chains.offsets[ci + 1];
    for (size_t k = cbi; k + 1 < cei; k++) {
      point p1 = hot_pixels[chains.indices[k]];
      point p2 = hot_pixels[chains.indices[k + 1]];
      auto x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
      auto x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
      auto seg_box = box{point{std::min(x1, x2), std::min(y1, y2)},
                         point{std::max(x1, x2), std::max(y1, y2)}};
      if (y1 < y2) {
        seg_data.push_back({x1, y1, x2, y2, 2 * ci});
        seg_boxes.push_back(seg_box);
      } else if (y2 < y1) {
        seg_data.push_back({x2, y2, x1, y1, 2 * ci + 1});
        seg_boxes.push_back(seg_box);
      }
    }
  }
  best_clipper::uniform_grid::grid seg_grid(seg_boxes);

  std::vector<std::size_t> exterior_half_chains;
  std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;

  for (std::size_t component_id = 0; component_id < num_components;
       ++component_id) {
    std::size_t leftmost_vertex = leftmost_vertexes[component_id];
    std::size_t start_hc{~0ULL};
    if (sorted_half_chains_offsets[leftmost_vertex] <
        sorted_half_chains_offsets[leftmost_vertex + 1]) {
      start_hc =
          sorted_half_chains[sorted_half_chains_offsets[leftmost_vertex]].id;
    } else {
      assert(false); // to do
    }

    auto min_x = bg::get<0>(hot_pixels[leftmost_vertex]);
    int32_t ray_y = bg::get<1>(hot_pixels[leftmost_vertex]);
    auto hit = cast_ray_minus_x(min_x, ray_y, seg_grid, seg_data);
    if (hit == (std::size_t)-1) {
      exterior_half_chains.push_back(start_hc);
    } else {
      ray_pairs.emplace_back(start_hc, hit);
    }
  }

  return fe_tuple{std::move(exterior_half_chains), std::move(ray_pairs)};
}

// ---------------------------------------------------------------------------
// compute_winding
// ---------------------------------------------------------------------------

std::vector<int> compute_winding(
    const chain_build_result &chains,
    const std::vector<std::pair<std::size_t, std::size_t>> &coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs,
    const std::vector<std::size_t> &exterior_half_chains) {

  auto t0 = std::chrono::high_resolution_clock::now();
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;

  std::vector<edge_with_power_t> edges;
  edges.reserve((coplanar_pairs.size() + ray_pairs.size()) * 2);
  for (auto [a, b] : coplanar_pairs) {
    edges.emplace_back(a, b, 0);
    edges.emplace_back(b, a, 0);
  }
  for (auto [a, b] : ray_pairs) {
    edges.emplace_back(a, b, 0);
    edges.emplace_back(b, a, 0);
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  auto [adj, adjacency] = bucket_sort(
      edges, num_half_chains,
      [](const edge_with_power_t &e) { return e.start; },
      [](const edge_with_power_t &e) { return std::pair{e.end, e.power}; });
  auto t2 = std::chrono::high_resolution_clock::now();

  constexpr int UNKNOWN = std::numeric_limits<int>::max() / 2;
  std::vector<int> winding(num_half_chains, UNKNOWN);
  for (auto exterior : exterior_half_chains)
    winding[exterior] = 0;

  std::vector<std::size_t> stack(exterior_half_chains.begin(),
                                 exterior_half_chains.end());
  while (!stack.empty()) {
    auto u = stack.back();
    stack.pop_back();
    for (auto j = adj[u]; j < adj[u + 1]; j++) {
      auto [v, diff] = adjacency[j];
      if (winding[v] == UNKNOWN) {
        winding[v] = winding[u] + diff;
        stack.push_back(v);
      }
    }
    // Dual edge: i ↔ i^1 with winding difference = -power of i
    auto dual = u ^ 1;
    if (winding[dual] == UNKNOWN) {
      winding[dual] = winding[u] - half_chain{u}.power(chains);
      stack.push_back(dual);
    }
  }
  auto t3 = std::chrono::high_resolution_clock::now();
  auto ms = [](auto d) {
    return std::chrono::duration<double, std::milli>(d).count();
  };
  std::fprintf(stderr,
               "  [winding] build=%.1fms sort=%.1fms dfs=%.1fms (hc=%zu)\n",
               ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), num_half_chains);
  return winding;
}

} // namespace best_clipper
