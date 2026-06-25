#include "filter_half_chains.hpp"
#include "include/graph_helper.hpp"

namespace best_clipper {

std::vector<int> compute_winding(
    const chain_group &chains,
    const half_chain_relations_t& half_chain_relations) {
  const auto& next_half_chain = half_chain_relations.next_half_chain;
  const auto& exterior_half_chains = half_chain_relations.exterior_half_chains;
  const auto& ray_pairs = half_chain_relations.ray_pairs;

  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;

  std::vector<edge_t> edges;
  edges.reserve(ray_pairs.size() * 2);
  for (auto [a, b] : ray_pairs) {
    edges.emplace_back(a.id, b.id);
    edges.emplace_back(b.id, a.id);
  }

  auto [adj, adjacency] = bucket_sort(
      edges, num_half_chains,
      [](const edge_t &e) { return e.start; },
      [](const edge_t &e) { return e.end; }
    );

  constexpr int UNKNOWN = std::numeric_limits<int>::max() / 2;
  std::vector<int> winding(num_half_chains, UNKNOWN);
  for (auto exterior : exterior_half_chains)
    winding[exterior.id] = 0;

  std::vector<half_chain> stack(exterior_half_chains.begin(),
                                 exterior_half_chains.end());
  while (!stack.empty()) {
    auto u = stack.back();
    stack.pop_back();
    for (auto j = adj[u.id]; j < adj[u.id + 1]; j++) {
      auto v = adjacency[j];
      if (winding[v] == UNKNOWN) {
        winding[v] = winding[u.id];
        stack.push_back(half_chain{v});
      }
    }
    // Dual edge: i ↔ i^1 with winding difference = -power of i
    auto dual = u.dual();
    if (winding[dual.id] == UNKNOWN) {
      winding[dual.id] = winding[u.id] - u.power(chains);
      stack.push_back(dual);
    }
    auto next = next_half_chain[u.id];
    if (winding[next.id] == UNKNOWN) {
      winding[next.id] = winding[u.id];
      stack.push_back(next);
    }
  }
  return winding;
}

}