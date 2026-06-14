#pragma once
#include "chain_builder.hpp"
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

struct half_chain {
  std::size_t id;
  std::size_t chain_id() const { return id / 2; }
  bool is_forward() const { return id % 2 == 0; }
  std::size_t source_node(const chain_build_result &c) const {
    if (is_forward())
      return c.indices[c.offsets[chain_id()]];
    else
      return c.indices[c.offsets[chain_id() + 1] - 1];
  }
  std::size_t target_node(const chain_build_result &c) const {
    if (is_forward())
      return c.indices[c.offsets[chain_id() + 1] - 1];
    else
      return c.indices[c.offsets[chain_id()]];
  }
  std::size_t next_along_source(const chain_build_result &c) const {
    if (is_forward())
      return c.indices[c.offsets[chain_id()] + 1];
    else
      return c.indices[c.offsets[chain_id() + 1] - 2];
  }
  int power(const chain_build_result &c) const {
    if (is_forward())
      return c.powers[chain_id()];
    else
      return -c.powers[chain_id()];
  }
  half_chain dual() const { return {id ^ 1}; }
};

using hcg_tuple = std::tuple<std::vector<half_chain>, std::vector<std::size_t>,
                             std::vector<half_chain>,
                             std::vector<std::pair<std::size_t, std::size_t>>>;

hcg_tuple build_half_chain_graph(const chain_build_result &chains,
                                 const std::vector<point> &hot_pixels);

} // namespace best_clipper
