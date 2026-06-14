#pragma once
#include "chain_builder.hpp"

#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

// Returns (offsets, data) where offsets has bucket_size+1 entries.
// Bucket i spans data[offsets[i] .. offsets[i+1]).
inline auto bucket_sort(auto vec, auto bucket_size, auto get_bucket,
                        auto get_left) {
  std::vector<std::size_t> times(bucket_size);
  for (auto val : vec)
    times[get_bucket(val)]++;
  std::vector<std::size_t> offsets(bucket_size + 1);
  std::exclusive_scan(std::begin(times), std::end(times), std::begin(offsets),
                      std::size_t{0});
  offsets[bucket_size] = vec.size();
  std::vector<std::invoke_result_t<decltype(get_left),
                                   typename decltype(vec)::value_type>>
      data(vec.size());
  auto cursors = offsets;
  for (auto val : vec)
    data[cursors[get_bucket(val)]++] = get_left(val);
  return std::pair{std::move(offsets), std::move(data)};
}

std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>> &edges);

chain_build_result
build_chains(const std::vector<edge_with_power_t> &sorted_edges,
             std::size_t node_num);

} // namespace best_clipper