#include "half_chain_graph.hpp"

#include <utility>
#include <vector>

namespace best_clipper {

half_chain_group_t build_half_chain_graph(const chain_group &chains) {
  assert(chains.out_offsets.size() == chains.in_offsets.size());
  std::size_t num_end_points = chains.out_offsets.size() - 1;
  std::size_t num_chains = chains.offsets.size() - 1;

  std::vector<std::size_t> bucket_half_chains_offsets(num_end_points + 1);
  for (std::size_t i = 0; i < num_end_points; ++i) {
    auto out_cnt = chains.out_offsets[i + 1] - chains.out_offsets[i];
    auto in_cnt = chains.in_offsets[i + 1] - chains.in_offsets[i];
    bucket_half_chains_offsets[i + 1] =
        bucket_half_chains_offsets[i] + out_cnt + in_cnt;
  }

  std::vector<half_chain_t> bucket_half_chains(num_chains * 2);
  auto cursors = bucket_half_chains_offsets;
  for (std::size_t i = 0; i < num_end_points; ++i) {
    for (auto ci = chains.out_offsets[i]; ci < chains.out_offsets[i + 1]; ++ci)
      bucket_half_chains[cursors[i]++] =
          half_chain_t{chains.out_chains[ci] * 2};
    for (auto ci = chains.in_offsets[i]; ci < chains.in_offsets[i + 1]; ++ci)
      bucket_half_chains[cursors[i]++] =
          half_chain_t{chains.in_chains[ci] * 2 + 1};
  }

  return {std::move(bucket_half_chains), std::move(bucket_half_chains_offsets)};
}

} // namespace best_clipper
