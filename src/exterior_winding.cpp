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
  int64_t x1, y1, dx, dy;
  std::size_t fwd_hc, rev_hc;
};

// ---------------------------------------------------------------------------
// cast_ray_minus_x
// ---------------------------------------------------------------------------

std::size_t cast_ray_minus_x(int32_t vx, int64_t ray_y,
                             const best_clipper::uniform_grid::grid &seg_grid,
                             const std::vector<chain_seg> &seg_data) {
  int64_t best_x = std::numeric_limits<int64_t>::min();
  std::size_t hit_id = (std::size_t)-1;
  int64_t hit_dy = 0;
  int64_t best_dx = 0;
  int64_t hit_y_low = 0;

  auto try_edge = [&](int64_t ix, int64_t dx, int64_t dy, int64_t y_low,
                      std::size_t half_chain_id) {
    if (ix >= vx)
      return;
    bool better = false;
    if (ix > best_x) {
      better = true;
    } else if (ix == best_x) {
      int sign_a = (dx ^ hit_dy) >= 0 ? 1 : -1;
      int sign_b = (best_dx ^ dy) >= 0 ? 1 : -1;
      bool same_sign = sign_a == sign_b;
      uint64_t abs_a = (uint64_t)(dx >= 0 ? dx : -dx) *
                       (uint64_t)(hit_dy >= 0 ? hit_dy : -hit_dy);
      uint64_t abs_b = (uint64_t)(best_dx >= 0 ? best_dx : -best_dx) *
                       (uint64_t)(dy >= 0 ? dy : -dy);
      if (same_sign ? (sign_a > 0 ? abs_a > abs_b : abs_a < abs_b)
                    : sign_a > sign_b) {
        better = true;
      } else if (same_sign && abs_a == abs_b) {
        uint64_t da = ray_y - y_low;
        uint64_t db = ray_y - hit_y_low;
        uint64_t abs_dy = dy > 0 ? (uint64_t)dy : (uint64_t)(-dy);
        uint64_t abs_hdy = hit_dy > 0 ? (uint64_t)hit_dy : (uint64_t)(-hit_dy);
        if (da * abs_hdy < db * abs_dy) {
          better = true;
        } else if (da * abs_hdy == db * abs_dy && abs_dy > abs_hdy) {
          better = true;
        }
      }
    }
    if (better) {
      best_x = ix;
      best_dx = dx;
      hit_dy = dy;
      hit_id = half_chain_id;
      hit_y_low = y_low;
    }
  };

  auto intersect_x = [&](int64_t x1, int64_t y1, int64_t dx,
                         int64_t dy) -> int64_t {
    return (int64_t)(((int128_t)x1 * dy + (int128_t)(ray_y - y1) * dx) / dy);
  };

  auto query_box =
      box{point{INT32_MIN, (int32_t)ray_y}, point{vx, (int32_t)ray_y}};
  seg_grid.query_intersects(query_box, [&](std::size_t idx) {
    auto &seg = seg_data[idx];
    int64_t y2 = seg.y1 + seg.dy;
    if (seg.y1 <= ray_y && ray_y < y2) {
      int64_t ix = intersect_x(seg.x1, seg.y1, seg.dx, seg.dy);
      try_edge(ix, seg.dx, seg.dy, seg.y1, seg.fwd_hc);
    } else if (y2 <= ray_y && ray_y < seg.y1) {
      int64_t ix = intersect_x(seg.x1, seg.y1, seg.dx, seg.dy);
      try_edge(ix, -seg.dx, -seg.dy, y2, seg.rev_hc);
    }
  });

  if (hit_id != (std::size_t)-1)
    return (hit_dy > 0) ? hit_id : (hit_id ^ 1);
  return (std::size_t)-1;
}

// ---------------------------------------------------------------------------
// find_exterior
// ---------------------------------------------------------------------------

fe_tuple find_exterior(const chain_build_result &chains,
                       const std::vector<point> &hot_pixels,
                       const std::vector<half_chain> &sorted_half_chains,
                       const std::vector<std::size_t> &half_chains) {
  std::size_t num_vertices = hot_pixels.size();

  std::vector<std::pair<std::size_t, std::size_t>> vertex_edges;
  for (auto h : sorted_half_chains)
    vertex_edges.emplace_back(h.source_node(chains), h.target_node(chains));
  for (std::size_t c = 0; c + 1 < chains.offsets.size(); c++) {
    auto chain_begin_idx = chains.offsets[c],
         chain_end_idx = chains.offsets[c + 1];
    for (std::size_t k = chain_begin_idx; k + 1 < chain_end_idx; k++)
      vertex_edges.emplace_back(chains.indices[k], chains.indices[k + 1]);
  }

  auto t_cc0 = std::chrono::high_resolution_clock::now();
  auto component_id = connected_components(num_vertices, vertex_edges);
  auto t_cc1 = std::chrono::high_resolution_clock::now();

  std::size_t num_components = 0;
  for (auto c : component_id)
    num_components = std::max(num_components, c + 1);
  std::vector<std::vector<std::size_t>> vertex_components(num_components);
  for (std::size_t v = 0; v < num_vertices; v++)
    vertex_components[component_id[v]].push_back(v);
  auto t_group = std::chrono::high_resolution_clock::now();
  auto ms = [](auto d) {
    return std::chrono::duration<double, std::milli>(d).count();
  };
  std::fprintf(stderr,
               "  [exterior] cc=%.1fms group=%.1fms (n=%zu m=%zu comps=%zu)\n",
               ms(t_cc1 - t_cc0), ms(t_group - t_cc1), num_vertices,
               vertex_edges.size(), num_components);

  auto t_data0 = std::chrono::high_resolution_clock::now();
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
      int32_t x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
      int32_t x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
      seg_data.push_back({(int64_t)x1, (int64_t)y1, (int64_t)x2 - (int64_t)x1,
                          (int64_t)y2 - (int64_t)y1, 2 * ci, 2 * ci + 1});
      seg_boxes.push_back(box{point{std::min(x1, x2), std::min(y1, y2)},
                              point{std::max(x1, x2), std::max(y1, y2)}});
    }
  }
  best_clipper::uniform_grid::grid seg_grid(seg_boxes);
  auto t_data1 = std::chrono::high_resolution_clock::now();

  std::vector<std::size_t> exterior_half_chains;
  std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;

  for (auto &vertex_component : vertex_components) {
    std::size_t leftmost_vertex = ~0ULL, first_hc_vertex = ~0ULL;
    int32_t min_x = std::numeric_limits<int32_t>::max();
    for (auto vertex : vertex_component) {
      int32_t x = bg::get<0>(hot_pixels[vertex]);
      if (x < min_x) {
        min_x = x;
        leftmost_vertex = vertex;
      }
      if (first_hc_vertex == ~0ULL &&
          half_chains[vertex] < half_chains[vertex + 1])
        first_hc_vertex = vertex;
    }
    if (first_hc_vertex == ~0ULL)
      continue;
    auto first_hc = sorted_half_chains[half_chains[first_hc_vertex]].id;

    int32_t ray_y = bg::get<1>(hot_pixels[leftmost_vertex]);
    auto hit = cast_ray_minus_x(min_x, ray_y, seg_grid, seg_data);
    if (hit == (std::size_t)-1) {
      exterior_half_chains.push_back(first_hc);
    } else {
      ray_pairs.emplace_back(first_hc, hit);
    }
  }

  auto t_rays = std::chrono::high_resolution_clock::now();
  std::fprintf(stderr, "  [exterior] data=%.1fms grid=%.1fms rays=%.1fms\n",
               ms(t_data0 - t_group), ms(t_data1 - t_data0),
               ms(t_rays - t_data1));
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
