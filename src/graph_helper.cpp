#include "include/graph_helper.hpp"

#include <cassert>

namespace best_clipper {

std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>> &edges) {

  std::vector<std::size_t> parent(n);
  std::vector<std::size_t> rank(n, 0);
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
    if (ra == rb)
      continue;
    if (rank[ra] < rank[rb])
      parent[ra] = rb;
    else if (rank[ra] > rank[rb])
      parent[rb] = ra;
    else {
      parent[ra] = rb;
      rank[rb]++;
    }
  }

  std::vector<std::size_t> comp(n);
  std::vector<std::size_t> root_to_id(n, ~0ULL);
  std::vector<std::size_t> comp_sizes;
  std::size_t num_components = 0;
  for (std::size_t v = 0; v < n; v++) {
    auto r = find(v);
    if (root_to_id[r] == ~0ULL) {
      root_to_id[r] = num_components++;
      comp_sizes.push_back(0);
    }
    auto cid = root_to_id[r];
    comp[v] = cid;
    comp_sizes[cid]++;
  }

  return {std::move(comp), std::move(comp_sizes)};
}

// ---------------------------------------------------------------------------
// build_chains chain decomposition from directed edges with power
// ---------------------------------------------------------------------------

std::tuple<std::vector<std::size_t>, std::vector<std::size_t>, std::vector<int>,
           std::vector<std::size_t>, std::vector<std::size_t>,
           std::vector<std::size_t>, std::vector<std::size_t>>
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
  for (std::size_t i = 0; i < node_num; i++) {
    is_end[i] =
        !(out_deg[i] == 1 && in_deg[i] == 1 && out_power[i] == in_power[i]);
  }

  std::vector<bool> visited(node_num), edge_used(sorted_edges.size());
  std::vector<std::size_t> idx, off{0};
  std::vector<int> powers;

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
      cur = sorted_edges[nj].end;
    }
    idx.push_back(cur);
    off.push_back(idx.size());
  }

  std::size_t num_chains = off.size() - 1;
  std::vector<std::size_t> chain_out_deg(node_num), chain_in_deg(node_num);
  for (std::size_t i = 0; i + 1 < off.size(); ++i) {
    chain_out_deg[idx[off[i]]]++;
    chain_in_deg[idx[off[i + 1] - 1]]++;
  }

  std::vector<std::size_t> full_out_offsets(node_num + 1, 0),
      out_chains(num_chains);
  std::vector<std::size_t> full_in_offsets(node_num + 1, 0),
      in_chains(num_chains);
  for (std::size_t i = 0; i < node_num; ++i) {
    full_out_offsets[i + 1] = full_out_offsets[i] + chain_out_deg[i];
  }
  for (std::size_t i = 0; i < node_num; ++i) {
    full_in_offsets[i + 1] = full_in_offsets[i] + chain_in_deg[i];
  }
  {
    auto out_cursors = full_out_offsets;
    auto in_cursors = full_in_offsets;
    for (std::size_t i = 0; i + 1 < off.size(); ++i) {
      auto start_idx = idx[off[i]];
      auto end_idx = idx[off[i + 1] - 1];
      out_chains[out_cursors[start_idx]++] = i;
      in_chains[in_cursors[end_idx]++] = i;
    }
  }
  {
    std::vector<std::size_t> compact_out{0}, compact_in{0};
    for (std::size_t i = 0; i < node_num; ++i) {
      if (full_in_offsets[i] < full_in_offsets[i + 1]) {
        compact_out.push_back(full_out_offsets[i + 1]);
        compact_in.push_back(full_in_offsets[i + 1]);
      }
    }
    full_out_offsets = std::move(compact_out);
    full_in_offsets = std::move(compact_in);
  }

  return {std::move(idx),        std::move(off),
          std::move(powers),     std::move(full_out_offsets),
          std::move(out_chains), std::move(full_in_offsets),
          std::move(in_chains)};
}

} // namespace best_clipper
