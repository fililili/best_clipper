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
  std::size_t source_node(const chain_group &chains) const {
    if (is_forward())
      return chains.indices[chains.offsets[chain_id()]];
    else
      return chains.indices[chains.offsets[chain_id() + 1] - 1];
  }
  std::size_t target_node(const chain_group &chains) const {
    if (is_forward())
      return chains.indices[chains.offsets[chain_id() + 1] - 1];
    else
      return chains.indices[chains.offsets[chain_id()]];
  }
  std::size_t next_along_source(const chain_group &chains) const {
    if (is_forward())
      return chains.indices[chains.offsets[chain_id()] + 1];
    else
      return chains.indices[chains.offsets[chain_id() + 1] - 2];
  }
  int power(const chain_group &chains) const {
    if (is_forward())
      return chains.powers[chain_id()];
    else
      return -chains.powers[chain_id()];
  }
  auto operator<=> (const half_chain_t &other) const = default;
  half_chain_t dual() const { return {id ^ 1}; }
};

using half_chain_group_t = std::tuple<std::vector<half_chain_t>, std::vector<std::size_t>>;

half_chain_group_t build_half_chain_graph(const chain_group &chains);

} // namespace best_clipper
