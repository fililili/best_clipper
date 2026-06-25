#pragma once
#include "chain_builder.hpp"
#include "half_chain_graph.hpp"
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

using fe_tuple = std::tuple<std::vector<half_chain>,
                            std::vector<std::size_t>,
                            std::vector<std::pair<std::size_t, std::size_t>>>;

fe_tuple build_half_chain_relations(const chain_group &chains,
                       const std::vector<point> &hot_pixels,
                       std::vector<half_chain> bucket_half_chains,
                       const std::vector<std::size_t> &bucket_half_chains_offsets);

} // namespace best_clipper
