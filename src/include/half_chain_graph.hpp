#pragma once
#include "chain_builder.hpp"
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

struct half_chain_t {
  std::size_t id;
  std::size_t chain_id() const { return id / 2; }
  bool is_forward() const { return id % 2 == 0; }
  point source_node(const chain_group_t &chains) const {
    if (is_forward())
      return chains.points[chains.offsets[chain_id()]];
    else
      return chains.points[chains.offsets[chain_id() + 1] - 1];
  }
  point target_node(const chain_group_t &chains) const {
    if (is_forward())
      return chains.points[chains.offsets[chain_id() + 1] - 1];
    else
      return chains.points[chains.offsets[chain_id()]];
  }
  point next_along_source(const chain_group_t &chains) const {
    if (is_forward())
      return chains.points[chains.offsets[chain_id()] + 1];
    else
      return chains.points[chains.offsets[chain_id() + 1] - 2];
  }
  int power(const chain_group_t &chains) const {
    if (is_forward())
      return chains.powers[chain_id()];
    else
      return -chains.powers[chain_id()];
  }
  auto operator<=>(const half_chain_t &other) const = default;
  half_chain_t dual() const { return {id ^ 1}; }
};

using half_chain_group_t =
    std::tuple<std::vector<half_chain_t>, std::vector<std::size_t>>;

half_chain_group_t
build_half_chain_graph(const chain_group_t &chains,
                       const std::vector<std::size_t> &out_offsets,
                       const std::vector<std::size_t> &out_chains,
                       const std::vector<std::size_t> &in_offsets,
                       const std::vector<std::size_t> &in_chains);

} // namespace best_clipper
