#pragma once
#include "chain_builder.hpp"
#include "graph_helper.hpp"
#include "half_chain_graph.hpp"
#include "connect_half_chains.hpp"

#include <utility>
#include <vector>

namespace best_clipper {

multi_polygon build_output(
    const chain_group &chains, const std::vector<point> &hot_pixels,
    half_chain_relations_t half_chain_relations, std::vector<bool> survive);

} // namespace best_clipper
