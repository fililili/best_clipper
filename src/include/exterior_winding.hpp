#pragma once
#include "chain_builder.hpp"
#include "half_chain_graph.hpp"
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

using fe_tuple = std::tuple<std::vector<std::size_t>,
                            std::vector<std::pair<std::size_t, std::size_t>>>;

fe_tuple find_exterior(const chain_build_result &chains,
                       const std::vector<point> &hot_pixels,
                       const std::vector<half_chain> &sorted_half_chains,
                       const std::vector<std::size_t> &half_chains);

std::vector<int> compute_winding(
    const chain_build_result &chains,
    const std::vector<std::pair<std::size_t, std::size_t>> &coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs,
    const std::vector<std::size_t> &exterior_half_chains);

} // namespace best_clipper
