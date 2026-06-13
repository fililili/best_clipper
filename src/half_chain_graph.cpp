#include "include/half_chain_graph.hpp"
#include "include/graph_helper.hpp"
#include "include/snap_rounding_helper.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace best_clipper {

namespace bg = boost::geometry;

hcg_tuple build_half_chain_graph(const chain_build_result &chains,
                                 const std::vector<point> &hot_pixels) {
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
  std::size_t num_vertices = hot_pixels.size();

  std::vector<half_chain> all;
  all.reserve(num_half_chains);
  for (std::size_t i = 0; i < num_half_chains; i++)
    all.push_back({i});

  auto [begin_loc, end_loc, sorted_half_chains] = bucket_sort(
      all, num_vertices, [&](half_chain h) { return h.source_node(chains); },
      [](half_chain h) { return h; });

  for (std::size_t v = 0; v < num_vertices; v++) {
    auto vertex_begin = begin_loc[v], vertex_end = end_loc[v];
    if (vertex_end - vertex_begin < 2)
      continue;
    std::sort(sorted_half_chains.begin() + vertex_begin,
              sorted_half_chains.begin() + vertex_end,
              [&](half_chain a, half_chain b) {
                return less_by_direction(
                    hot_pixels[v], hot_pixels[a.next_along_source(chains)],
                    hot_pixels[b.next_along_source(chains)]);
              });
  }

  std::vector<half_chain> next_half_chain(num_half_chains, {~0ULL});
  std::vector<std::pair<std::size_t, std::size_t>> coplanar;

  for (std::size_t v = 0; v < num_vertices; v++) {
    auto vertex_begin = begin_loc[v], vertex_end = end_loc[v];
    if (vertex_begin == vertex_end)
      continue;
    assert(vertex_begin + 1 != vertex_end);
    for (auto it = vertex_begin + 1; it < vertex_end; ++it) {
      auto prev = sorted_half_chains[it - 1], cur = sorted_half_chains[it];
      next_half_chain[prev.dual().id] = cur;
      coplanar.emplace_back(cur.id, prev.dual().id);
    }
    auto first = sorted_half_chains[vertex_begin],
         last = sorted_half_chains[vertex_end - 1];
    next_half_chain[last.dual().id] = first;
    coplanar.emplace_back(first.id, last.dual().id);
  }

  return hcg_tuple{std::move(sorted_half_chains),
                   std::vector<std::size_t>(begin_loc.begin(), begin_loc.end()),
                   std::vector<std::size_t>(end_loc.begin(), end_loc.end()),
                   std::move(next_half_chain), std::move(coplanar)};
}

} // namespace best_clipper
