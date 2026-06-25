#include "include/half_chain_graph.hpp"
#include "include/graph_helper.hpp"
#include "include/snap_rounding_helper.hpp"

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

namespace best_clipper {

hcg_tuple build_half_chain_graph(const chain_group &chains,
                                 const std::vector<point> &hot_pixels) {
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
  std::size_t num_vertices = hot_pixels.size();

  std::vector<half_chain> all;
  all.reserve(num_half_chains);
  for (std::size_t i = 0; i < num_half_chains; i++)
    all.push_back({i});

  auto [bucket_half_chains_offsets, bucket_half_chains] = bucket_sort(
      all, num_vertices, [&](half_chain h) { return h.source_node(chains); },
      [](half_chain h) { return h; });

  return hcg_tuple{std::move(bucket_half_chains), std::move(bucket_half_chains_offsets)};
}

} // namespace best_clipper
