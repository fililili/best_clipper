#include "include/graph_helper.hpp"

#include <cassert>
#include <cstdint>

namespace best_clipper {

std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>> &edges) {

  std::vector<std::size_t> parent(n);
  for (std::size_t i = 0; i < n; i++)
    parent[i] = i;

  auto find = [&](std::size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };

  for (auto [a, b] : edges) {
    auto ra = find(a), rb = find(b);
    if (ra != rb)
      parent[ra] = rb;
  }

  std::vector<std::size_t> comp(n);
  std::vector<std::size_t> root_to_id(n, ~0ULL);
  std::size_t num_components = 0;
  for (std::size_t v = 0; v < n; v++) {
    auto r = find(v);
    if (root_to_id[r] == ~0ULL)
      root_to_id[r] = num_components++;
    comp[v] = root_to_id[r];
  }

  return comp;
}

// ---------------------------------------------------------------------------
// build_chains — chain decomposition from directed edges with power
// ---------------------------------------------------------------------------

chain_build_result
build_chains(const std::vector<edge_with_power_t> &sorted_edges,
             std::size_t node_num) {
  std::vector<std::size_t> edge_offsets(node_num + 1, 0);
  {
    std::vector<std::size_t> edge_count(node_num, 0);
    for (auto &e : sorted_edges)
      edge_count[e.start]++;
    std::size_t cur = 0;
    for (std::size_t v = 0; v < node_num; v++) {
      edge_offsets[v] = cur;
      cur += edge_count[v];
    }
    edge_offsets.back() = cur;
  }

  std::vector<std::uint32_t> out_deg(node_num), in_deg(node_num);
  std::vector<int> out_power(node_num), in_power(node_num);
  for (const auto &e : sorted_edges) {
    out_deg[e.start]++;
    in_deg[e.end]++;
    out_power[e.start] = e.power;
    in_power[e.end] = e.power;
  }

  std::vector<bool> is_end(node_num);
  for (std::size_t i = 0; i < node_num; i++)
    is_end[i] =
        !(out_deg[i] == 1 && in_deg[i] == 1 && out_power[i] == in_power[i]);

  std::vector<bool> visited(node_num), edge_used(sorted_edges.size());
  std::vector<std::size_t> idx, off{0};
  std::vector<int> powers;
  std::vector<std::size_t> edge_to_chain(sorted_edges.size(), ~0ULL);

  for (std::size_t i = 0; i < node_num; i++) {
    if (!is_end[i])
      continue;
    visited[i] = true;
    for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++) {
      if (edge_used[j] || sorted_edges[j].start != i)
        continue;
      edge_used[j] = true;
      idx.push_back(i);
      powers.push_back(sorted_edges[j].power);
      edge_to_chain[j] = off.size() - 1;
      auto cur = sorted_edges[j].end;
      while (!is_end[cur]) {
        visited[cur] = true;
        idx.push_back(cur);
        std::size_t nj = ~0ULL;
        for (std::size_t k = edge_offsets[cur]; k < edge_offsets[cur + 1]; k++)
          if (!edge_used[k] && sorted_edges[k].start == cur) {
            nj = k;
            break;
          }
        assert(nj != ~0ULL);
        edge_used[nj] = true;
        edge_to_chain[nj] = off.size() - 1;
        cur = sorted_edges[nj].end;
      }
      idx.push_back(cur);
      off.push_back(idx.size());
    }
  }

  for (std::size_t i = 0; i < node_num; i++) {
    if (visited[i])
      continue;
    std::size_t sj = ~0ULL;
    for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++)
      if (!edge_used[j] && sorted_edges[j].start == i) {
        sj = j;
        break;
      }
    if (sj == ~0ULL)
      continue;
    visited[i] = true;
    edge_used[sj] = true;
    idx.push_back(i);
    powers.push_back(sorted_edges[sj].power);
    edge_to_chain[sj] = off.size() - 1;
    auto cur = sorted_edges[sj].end;
    while (i != cur) {
      visited[cur] = true;
      idx.push_back(cur);
      std::size_t nj = ~0ULL;
      for (std::size_t k = edge_offsets[cur]; k < edge_offsets[cur + 1]; k++)
        if (!edge_used[k] && sorted_edges[k].start == cur) {
          nj = k;
          break;
        }
      assert(nj != ~0ULL);
      edge_used[nj] = true;
      edge_to_chain[nj] = off.size() - 1;
      cur = sorted_edges[nj].end;
    }
    idx.push_back(cur);
    off.push_back(idx.size());
  }

  return {std::move(idx), std::move(off), std::move(powers),
          std::move(edge_to_chain)};
}

} // namespace best_clipper
