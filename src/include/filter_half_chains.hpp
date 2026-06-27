#pragma once

#include "chain_builder.hpp"
#include "connect_half_chains.hpp"
#include "half_chain_graph.hpp"

namespace best_clipper {

std::vector<int>
compute_winding(const chain_group_t &chains,
                const half_chain_relations_t &half_chain_relations);

inline std::vector<bool> filter_survive(const std::vector<int> &winding,
                                        auto filter_fn) {
  std::vector<bool> survive(winding.size());
  for (std::size_t i = 0; i < winding.size(); i++)
    survive[i] = filter_fn(winding[i]);
  return survive;
}
} // namespace best_clipper