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

  seg_grid.query_ray_left(vx, ray_y, [&](std::size_t idx) {
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
    return std::min(best_hit_x1, best_hit_x2);
    // we want to return best_hit_x_floor that satisfied best_hit_x_real >= best_hit_x_floor, 
    // and best_hit_x_real >= min(best_hit_x1, best_hit_x2), so use min(best_hit_x1, best_hit_x2) as the return value
  });

  return hit_half_chain_id;
}
constexpr auto less_by_direction_neg_x_split = [](point source, point target1, point target2) {
    enum class quadrant { _3, _4, _1, _2, zero };
    constexpr auto get_quadrant = [](int32_t dx, int32_t dy) -> quadrant {
        if (dx < 0 && dy <= 0) return quadrant::_3;
        else if (dx >= 0 && dy < 0) return quadrant::_4;
        else if (dx > 0 && dy >= 0) return quadrant::_1;
        else if (dx <= 0 && dy > 0) return quadrant::_2;
        return quadrant::zero;
    };
    
    int32_t dx1 = bg::get<0>(target1) - bg::get<0>(source);
    int32_t dy1 = bg::get<1>(target1) - bg::get<1>(source);
    int32_t dx2 = bg::get<0>(target2) - bg::get<0>(source);
    int32_t dy2 = bg::get<1>(target2) - bg::get<1>(source);

    auto q1 = get_quadrant(dx1, dy1);
    auto q2 = get_quadrant(dx2, dy2);
        if (q1 != q2) return q1 < q2;

    // Same quadrant: compare slopes using cross product.
    // slope1 < slope2  ⟺  cross = dy1*dx2 - dy2*dx1 < 0 in all quadrants.
    int64_t cross = (int64_t)dy1 * dx2 - (int64_t)dy2 * dx1;
    return cross < 0;
};
// ---------------------------------------------------------------------------
// find_exterior
// ---------------------------------------------------------------------------

fe_tuple
find_exterior(const chain_group &chains,
              const std::vector<point> &hot_pixels,
              const std::vector<half_chain> &sorted_half_chains,
              const std::vector<std::size_t> &sorted_half_chains_offsets) {
  std::vector<std::pair<std::size_t, std::size_t>> chain_connected_edges;
  for (std::size_t i = 0; i < sorted_half_chains_offsets.size() - 1; ++i) {
    auto begin_idx = sorted_half_chains_offsets[i],
         end_idx = sorted_half_chains_offsets[i + 1];
    for (auto k = begin_idx + 1; k < end_idx; ++k) {
      chain_connected_edges.emplace_back(sorted_half_chains[begin_idx].chain_id(),
                                         sorted_half_chains[k].chain_id());
    }
  }
  auto chain_component_ids =
      connected_components(chains.offsets.size() - 1, chain_connected_edges);
  std::size_t num_components = 0;
  for (auto c : chain_component_ids) {
    num_components = std::max(num_components, c + 1);
  }
  // todo: if num_components == 1, no need to create seg_data, seg_grid. Just return one exterior_half_chain and empty ray_pairs

  std::vector<point> ray_start_points(num_components);
  std::vector<std::size_t> ray_start_half_chains(num_components);
  {
    std::vector<coordinate_type> min_xs(num_components,
                                        std::numeric_limits<int32_t>::max());
    std::vector<std::size_t> chain_ids(num_components);
    std::vector<std::size_t> position_in_chains(num_components);
    for (std::size_t chain_id = 0; chain_id < chains.offsets.size() - 1;
        ++chain_id) {
      auto component_id = chain_component_ids[chain_id];
      auto chain_begin_idx = chains.offsets[chain_id],
          chain_end_idx = chains.offsets[chain_id + 1];
      for (std::size_t k = chain_begin_idx; k < chain_end_idx - 1; k++) {
        auto vertex = chains.indices[k];
        auto x = bg::get<0>(hot_pixels[vertex]);
        if (x < min_xs[component_id]) {
          min_xs[component_id] = x;
          chain_ids[component_id] = chain_id;
          position_in_chains[component_id] = k;
        }
      }
    }
    for (std::size_t component_id = 0; component_id < num_components; ++component_id) {
      auto vertex = chains.indices[position_in_chains[component_id]];
      ray_start_points[component_id] = hot_pixels[vertex];
      if (sorted_half_chains_offsets[vertex] <
          sorted_half_chains_offsets[vertex + 1]) {
        ray_start_half_chains[component_id] =
            sorted_half_chains[sorted_half_chains_offsets[vertex]].id;
      } else {
        auto prev_vertex = chains.indices[position_in_chains[component_id] - 1];
        auto next_vertex = chains.indices[position_in_chains[component_id] + 1];
        if (less_by_direction_neg_x_split(hot_pixels[vertex], hot_pixels[prev_vertex], hot_pixels[next_vertex])) {
            ray_start_half_chains[component_id] = 2 * chain_ids[component_id] + 1;
        } else {
            ray_start_half_chains[component_id] = 2 * chain_ids[component_id];
        }
      }
    }
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
    auto min_x = bg::get<0>(ray_start_points[component_id]);
    int32_t ray_y = bg::get<1>(ray_start_points[component_id]);
    auto hit = cast_ray_minus_x(min_x, ray_y, seg_grid, seg_data);
    if (hit == (std::size_t)-1) {
      exterior_half_chains.push_back(ray_start_half_chains[component_id]);
    } else {
      ray_pairs.emplace_back(ray_start_half_chains[component_id], hit);
    }
  }

  return fe_tuple{std::move(exterior_half_chains), std::move(ray_pairs)};
}

// ---------------------------------------------------------------------------
// compute_winding
// ---------------------------------------------------------------------------

std::vector<int> compute_winding(
    const chain_group &chains,
    const std::vector<std::pair<std::size_t, std::size_t>> &coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs,
    const std::vector<std::size_t> &exterior_half_chains) {
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

  auto [adj, adjacency] = bucket_sort(
      edges, num_half_chains,
      [](const edge_with_power_t &e) { return e.start; },
      [](const edge_with_power_t &e) { return std::pair{e.end, e.power}; });

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
    auto dual = half_chain{u}.dual().id;
    if (winding[dual] == UNKNOWN) {
      winding[dual] = winding[u] - half_chain{u}.power(chains);
      stack.push_back(dual);
    }
  }
  return winding;
}

} // namespace best_clipper
