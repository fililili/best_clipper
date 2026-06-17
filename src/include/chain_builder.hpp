#pragma once
#include "geometry_types.hpp"
#include "graph_types.hpp"
#include <tuple>
#include <vector>

namespace best_clipper {

struct chain_group {
  // todo: store coordinates direclty
  std::vector<std::size_t> indices, offsets;
  std::vector<int> powers;
};

std::tuple<std::vector<point>, chain_group>
build_chains_from_input(const std::vector<point> &points,
                        const std::vector<std::size_t> &offsets);

} // namespace best_clipper
