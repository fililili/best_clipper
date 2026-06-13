#pragma once
#include "chain_builder.hpp"

#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace best_clipper {

inline auto bucket_sort(auto vec, auto bucket_size, auto get_bucket,
                        auto get_left) {
  std::vector<std::size_t> times(bucket_size);
  for (auto val : vec)
    times[get_bucket(val)]++;
  std::vector<std::size_t> begin_location(times.size());
  std::exclusive_scan(std::begin(times), std::end(times),
                      std::begin(begin_location), std::size_t{0});
  std::vector<std::invoke_result_t<decltype(get_left),
                                   typename decltype(vec)::value_type>>
      left(vec.size());
  auto current_location{begin_location};
  for (auto val : vec)
    left[current_location[get_bucket(val)]++] = get_left(val);
  return std::tuple{std::move(begin_location), std::move(current_location),
                    std::move(left)};
}

std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>> &edges);

chain_build_result
build_chains(const std::vector<edge_with_power_t> &sorted_edges,
             std::size_t node_num);

} // namespace best_clipper