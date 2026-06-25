#include "include/output_builder.hpp"
#include "include/geometry_types.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <limits>

namespace best_clipper {

multi_polygon build_output(
    const chain_group &chains, const std::vector<point> &hot_pixels,
    half_chain_relations_t half_chain_relations, std::vector<bool> survive) {
  auto& next_half_chain = half_chain_relations.next_half_chain;
  const auto& exterior_half_chains = half_chain_relations.exterior_half_chains;
  const auto& ray_pairs = half_chain_relations.ray_pairs;

  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
  std::size_t num_chains = chains.offsets.size() - 1;

  // Build connected components: next_half_chain + ray + dual cancellation → same face
  std::vector<std::pair<std::size_t, std::size_t>> same_face_half_chains;
  same_face_half_chains.reserve(next_half_chain.size() + ray_pairs.size() + num_chains);

  for(std::size_t i = 0; i < next_half_chain.size(); i++) {
    same_face_half_chains.emplace_back(i, next_half_chain[i].id);
  }
  for(auto [start, end] : ray_pairs) {
    same_face_half_chains.emplace_back(start.id, end.id);
  }

  for (std::size_t chain_idx = 0; chain_idx < num_chains; chain_idx++) {
    std::size_t forward_id = 2 * chain_idx, reverse_id = 2 * chain_idx + 1;
    if (survive[forward_id] && survive[reverse_id]) {
      survive[forward_id] = survive[reverse_id] = false;
      same_face_half_chains.emplace_back(forward_id, reverse_id);
    }
  }

  for (std::size_t i = 0; i < num_half_chains; i++) {
    if (!survive[i])
      continue;
    while (!survive[next_half_chain[i].id]) {
      next_half_chain[i] = next_half_chain[next_half_chain[i].dual().id];
    }
  }

  auto component_id = connected_components(num_half_chains, same_face_half_chains);

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
  return result;
}

} // namespace best_clipper
