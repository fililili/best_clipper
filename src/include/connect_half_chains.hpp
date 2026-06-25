#pragma once
#include "chain_builder.hpp"
#include "half_chain_graph.hpp"
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

struct half_chain_relations_t {
    struct ray_pair {
        half_chain_t ray_start;
        half_chain_t ray_hit;
    };
    std::vector<half_chain_t> next_half_chain;
    std::vector<half_chain_t> exterior_half_chains;
    std::vector<ray_pair> ray_pairs;
};

half_chain_relations_t build_half_chain_relations(const chain_group &chains, std::vector<half_chain_t> bucket_half_chains,
                       const std::vector<std::size_t> &bucket_half_chains_offsets);

} // namespace best_clipper
