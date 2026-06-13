#pragma once
#include "chain_builder.hpp"
#include "graph_helper.hpp"
#include "half_chain_graph.hpp"

#include <utility>
#include <vector>

namespace best_clipper {

multi_polygon build_output(
    const chain_build_result &chains, const std::vector<point> &hot_pixels,
    std::vector<half_chain> next_half_chain, std::vector<bool> survive,
    const std::vector<std::pair<std::size_t, std::size_t>> &coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>> &ray_pairs);

} // namespace best_clipper
