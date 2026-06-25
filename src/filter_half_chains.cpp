#include "filter_half_chains.hpp"
#include "include/graph_helper.hpp"

namespace best_clipper {

std::vector<int> compute_winding(
    const chain_group &chains,
    const std::vector<half_chain> &next_half_chain,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs,
    const std::vector<std::size_t> &exterior_half_chains) {
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;

  std::vector<edge_t> edges;
  edges.reserve(ray_pairs.size() * 2);
  for (auto [a, b] : ray_pairs) {
    edges.emplace_back(a, b);
    edges.emplace_back(b, a);
  }

  auto [adj, adjacency] = bucket_sort(
      edges, num_half_chains,
      [](const edge_t &e) { return e.start; },
      [](const edge_t &e) { return e.end; }
    );

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
      auto v = adjacency[j];
      if (winding[v] == UNKNOWN) {
        winding[v] = winding[u];
        stack.push_back(v);
      }
    }
    // Dual edge: i ↔ i^1 with winding difference = -power of i
    auto dual = half_chain{u}.dual().id;
    if (winding[dual] == UNKNOWN) {
      winding[dual] = winding[u] - half_chain{u}.power(chains);
      stack.push_back(dual);
    }
    auto next = next_half_chain[u].id;
    if (winding[next] == UNKNOWN) {
      winding[next] = winding[u];
      stack.push_back(next);
    }
  }
  return winding;
}

}