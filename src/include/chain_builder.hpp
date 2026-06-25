#pragma once
#include "geometry_types.hpp"
#include "graph_types.hpp"
#include <tuple>
#include <vector>

namespace best_clipper {

struct chain_group {
  std::vector<point> points;
  std::vector<std::size_t> offsets;
  std::vector<int> powers;
};

std::tuple<chain_group, std::vector<std::size_t>, std::vector<std::size_t>,
           std::vector<std::size_t>, std::vector<std::size_t>>
build_chains_from_input(const std::vector<point> &points,
                        const std::vector<std::size_t> &offsets);

} // namespace best_clipper
