#pragma once

#include "chain_builder.hpp"
#include "half_chain_graph.hpp"

namespace best_clipper {

std::vector<int> compute_winding(
    const chain_group &chains,
    const std::vector<half_chain> &next_half_chain,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs,
    const std::vector<std::size_t> &exterior_half_chains);
    
inline std::vector<bool> filter_survive(const std::vector<int> &winding,
                                        auto filter_fn) {
  std::vector<bool> survive(winding.size());
  for (std::size_t i = 0; i < winding.size(); i++)
    survive[i] = filter_fn(winding[i]);
  return survive;
}
}