#pragma once
#include "geometry_types.hpp"
#include "graph_types.hpp"
#include <tuple>
#include <vector>

namespace best_clipper {

struct chain_build_result {
  std::vector<std::size_t> indices, offsets;
  std::vector<int> powers;
};

std::tuple<std::vector<point>, chain_build_result>
build_chains_from_input(const std::vector<point> &points,
                        const std::vector<std::size_t> &offsets);

} // namespace best_clipper
