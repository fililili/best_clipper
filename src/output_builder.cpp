#include "include/output_builder.hpp"
#include "include/geometry_types.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <limits>

namespace best_clipper {

multi_polygon build_output(
    const chain_build_result &chains, const std::vector<point> &hot_pixels,
    std::vector<half_chain> next_half_chain, std::vector<bool> survive,
    const std::vector<std::pair<std::size_t, std::size_t>> &coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs) {
  using clock = std::chrono::high_resolution_clock;
  auto t0 = clock::now();
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
  std::size_t num_chains = chains.offsets.size() - 1;

  // Build connected components: coplanar + ray + dual cancellation → same face
  std::vector<std::pair<std::size_t, std::size_t>> face_edges;
  face_edges.reserve(coplanar_pairs.size() + ray_pairs.size() + num_chains);
  face_edges.insert(face_edges.end(), coplanar_pairs.begin(),
                    coplanar_pairs.end());
  face_edges.insert(face_edges.end(), ray_pairs.begin(), ray_pairs.end());

  for (std::size_t chain_idx = 0; chain_idx < num_chains; chain_idx++) {
    std::size_t forward_id = 2 * chain_idx, reverse_id = 2 * chain_idx + 1;
    if (survive[forward_id] && survive[reverse_id]) {
      survive[forward_id] = survive[reverse_id] = false;
      face_edges.emplace_back(forward_id, reverse_id);
    }
  }

  for (std::size_t i = 0; i < num_half_chains; i++) {
    if (!survive[i])
      continue;
    while (!survive[next_half_chain[i].id]) {
      next_half_chain[i] = next_half_chain[next_half_chain[i].dual().id];
    }
  }
  auto t1 = clock::now();

  auto component_id = connected_components(num_half_chains, face_edges);
  auto t2 = clock::now();

  std::vector<std::pair<std::size_t, std::size_t>> half_chain_to_face;
  for (std::size_t i = 0; i < num_half_chains; i++) {
    if (survive[i]) {
      half_chain_to_face.emplace_back(component_id[i], i);
    }
  }

  std::size_t num_faces = 0;
  for (auto c : component_id)
    num_faces = std::max(num_faces, c + 1);
  auto [face_chains, existed_half_chains] = bucket_sort(
      half_chain_to_face, num_faces,
      [](const std::pair<std::size_t, std::size_t> &p) { return p.first; },
      [](const std::pair<std::size_t, std::size_t> &p) { return p.second; });
  auto t3 = clock::now();

  // For each face, trace its half-chains into a polygon
  multi_polygon result;
  std::vector<bool> done(num_half_chains, false);
  for (std::size_t f = 0; f < num_faces; f++) {
    std::vector<ring> rings;

    for (std::size_t j = face_chains[f]; j < face_chains[f + 1]; j++) {
      ring current_ring;
      std::size_t current_id = existed_half_chains[j];
      if (done[current_id])
        continue;
      do {
        auto h = half_chain{current_id};
        auto chain_idx = h.chain_id();
        auto chain_begin = chains.offsets[chain_idx],
             chain_end = chains.offsets[chain_idx + 1];
        if (h.is_forward()) {
          for (std::size_t k = chain_begin; k < chain_end - 1; k++)
            current_ring.push_back(hot_pixels[chains.indices[k]]);
        } else {
          for (std::size_t k = chain_end - 1; k > chain_begin; k--)
            current_ring.push_back(hot_pixels[chains.indices[k]]);
        }
        done[current_id] = true;
        current_id = next_half_chain[current_id].id;
      } while (current_id != existed_half_chains[j]);

      if (!current_ring.empty())
        current_ring.push_back(current_ring.front());
      rings.push_back(std::move(current_ring));
    }

    if (rings.empty())
      continue;

    std::size_t outer_index = 0;
    {
      int32_t min_x = std::numeric_limits<int32_t>::max();
      for (std::size_t ring_index = 0; ring_index < rings.size();
           ring_index++) {
        for (const auto &point_value : rings[ring_index]) {
          int32_t x = bg::get<0>(point_value);
          if (x < min_x) {
            min_x = x;
            outer_index = ring_index;
          }
        }
      }
    }

    polygon polygon_result;
    polygon_result.outer() = std::move(rings[outer_index]);
    for (std::size_t ring_index = 0; ring_index < rings.size(); ring_index++) {
      if (ring_index == outer_index || rings[ring_index].empty())
        continue;
      polygon_result.inners().push_back(std::move(rings[ring_index]));
    }

    result.push_back(std::move(polygon_result));
  }

  auto t4 = clock::now();
  auto ms = [](auto d) {
    return std::chrono::duration<double, std::milli>(d).count();
  };
  std::fprintf(stderr,
               "  [output] cancel=%.1fms cc=%.1fms bucket=%.1fms trace=%.1fms "
               "(faces=%zu)\n",
               ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), ms(t4 - t3), num_faces);
  return result;
}

} // namespace best_clipper
