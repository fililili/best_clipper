#include "include/half_chain_graph.hpp"
#include "include/graph_helper.hpp"
#include "include/snap_rounding_helper.hpp"

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

namespace best_clipper {

namespace bg = boost::geometry;

void sort_bucket_half_chains(const chain_group &chains,
                        const std::vector<point> &hot_pixels,
                        std::vector<half_chain> &bucket_half_chains,
                        std::vector<std::size_t> &bucket_half_chains_offsets) {
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
  std::size_t num_vertices = hot_pixels.size();

  for (std::size_t v = 0; v < num_vertices; v++) {
    auto vertex_begin = bucket_half_chains_offsets[v], vertex_end = bucket_half_chains_offsets[v + 1];
    auto n = vertex_end - vertex_begin;
    if (n < 2)
      continue;
    auto begin_it = bucket_half_chains.begin() + vertex_begin;

    point vpt = hot_pixels[v];
    int32_t vx = bg::get<0>(vpt), vy = bg::get<1>(vpt);

    // Bucket by octant (8 buckets, octant 0 starts at -x).
    std::array<std::size_t, 8> bucket_cnt{};
    std::vector<int> octants(n);
    for (std::size_t i = 0; i < n; i++) {
      auto h = begin_it[i];
      auto pt = hot_pixels[h.next_along_source(chains)];
      int64_t dx = (int64_t)bg::get<0>(pt) - vx;
      int64_t dy = (int64_t)bg::get<1>(pt) - vy;
      int oct;
      if (dx <= 0) {
        if (dy <= 0)
          oct = (-dx >= -dy) ? 0 : 1;
        else
          oct = (-dx >= dy) ? 7 : 6;
      } else {
        if (dy <= 0)
          oct = (dx >= -dy) ? 3 : 2;
        else
          oct = (dx >= dy) ? 4 : 5;
      }
      octants[i] = oct;
      bucket_cnt[oct]++;
    }

    // Compute bucket start positions and scatter.
    std::array<std::size_t, 9> bucket_pos{};
    for (int o = 0; o < 8; o++)
      bucket_pos[o + 1] = bucket_pos[o] + bucket_cnt[o];
    std::vector<half_chain> temp(n);
    auto cur_pos = bucket_pos;
    for (std::size_t i = 0; i < n; i++)
      temp[cur_pos[octants[i]]++] = begin_it[i];

    // Sort each bucket by cross product.
    for (int o = 0; o < 8; o++) {
      std::sort(temp.begin() + bucket_pos[o], temp.begin() + bucket_pos[o + 1],
                [&](half_chain a, half_chain b) {
                  auto a_pt = hot_pixels[a.next_along_source(chains)];
                  auto b_pt = hot_pixels[b.next_along_source(chains)];
                  int64_t ax = (int64_t)bg::get<0>(a_pt) - vx;
                  int64_t ay = (int64_t)bg::get<1>(a_pt) - vy;
                  int64_t bx = (int64_t)bg::get<0>(b_pt) - vx;
                  int64_t by = (int64_t)bg::get<1>(b_pt) - vy;
                  return ax * by - ay * bx > 0;
                });
    }

    std::move(temp.begin(), temp.end(), begin_it);
  }
}

std::vector<half_chain> build_next_chains(std::vector<half_chain> &bucket_half_chains, std::vector<std::size_t> &bucket_half_chains_offsets) {
  std::vector<half_chain> next_half_chain(bucket_half_chains.size(), {~0ULL});

  for (std::size_t i = 0; i + 1 < bucket_half_chains_offsets.size(); ++i) {
    auto vertex_begin = bucket_half_chains_offsets[i], vertex_end = bucket_half_chains_offsets[i + 1];
    if (vertex_begin == vertex_end) {
      continue; // skip by hand for head tail next relation
    }
    for (auto it = vertex_begin + 1; it < vertex_end; ++it) {
      auto prev = bucket_half_chains[it - 1], cur = bucket_half_chains[it];
      next_half_chain[prev.dual().id] = cur;
    }
    auto first = bucket_half_chains[vertex_begin],
         last = bucket_half_chains[vertex_end - 1];
    next_half_chain[last.dual().id] = first;
  }
  return next_half_chain;
}

hcg_tuple build_half_chain_graph(const chain_group &chains,
                                 const std::vector<point> &hot_pixels) {
  std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
  std::size_t num_vertices = hot_pixels.size();

  std::vector<half_chain> all;
  all.reserve(num_half_chains);
  for (std::size_t i = 0; i < num_half_chains; i++)
    all.push_back({i});

  auto [bucket_half_chains_offsets, bucket_half_chains] = bucket_sort(
      all, num_vertices, [&](half_chain h) { return h.source_node(chains); },
      [](half_chain h) { return h; });
      
  sort_bucket_half_chains(chains, hot_pixels, bucket_half_chains, bucket_half_chains_offsets);

  auto next_half_chain = build_next_chains(bucket_half_chains, bucket_half_chains_offsets);

  return hcg_tuple{std::move(bucket_half_chains), std::move(bucket_half_chains_offsets),
                   std::move(next_half_chain)};
}

} // namespace best_clipper
