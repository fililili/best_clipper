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
// build_chains — chain decomposition from directed edges with power.
// Merges adjacent edges that share a vertex with out_deg==in_deg==1 and same
// power. Returns (chain_idx, chain_off, chain_powers,
//                out_offsets, out_chains, in_offsets, in_chains).
// ---------------------------------------------------------------------------

std::tuple<std::vector<std::size_t>, std::vector<std::size_t>, std::vector<int>,
           std::vector<std::size_t>, std::vector<std::size_t>,
           std::vector<std::size_t>, std::vector<std::size_t>>
build_chains(const std::vector<edge_with_power_t> &sorted_edges,
             std::size_t node_num) {
  // 1. Degrees and powers.
  std::vector<std::uint32_t> out_deg(node_num), in_deg(node_num);
  std::vector<int> out_power(node_num), in_power(node_num);
  for (const auto &e : sorted_edges) {
    out_deg[e.start]++;
    in_deg[e.end]++;
    out_power[e.start] = e.power;
    in_power[e.end] = e.power;
  }

  // 2. Edge offsets from out_deg.
  std::vector<std::size_t> edge_offsets(node_num + 1, 0);
  std::exclusive_scan(out_deg.begin(), out_deg.end(), edge_offsets.begin(),
                      std::size_t{0});
  edge_offsets.back() = sorted_edges.size();

  // 3. End-node classification.
  std::vector<bool> is_end(node_num);
  for (std::size_t i = 0; i < node_num; i++)
    is_end[i] =
        !(out_deg[i] == 1 && in_deg[i] == 1 && out_power[i] == in_power[i]);

  // 4. Compact mapping and offsets from end nodes.
  std::vector<std::size_t> compact_of(node_num, ~0ULL);
  std::vector<std::size_t> out_offsets{0}, in_offsets{0};
  for (std::size_t v = 0; v < node_num; v++) {
    if (!is_end[v] || in_deg[v] == 0)
      continue;
    compact_of[v] = out_offsets.size() - 1;
    out_offsets.push_back(out_offsets.back() + (std::size_t)out_deg[v]);
    in_offsets.push_back(in_offsets.back() + (std::size_t)in_deg[v]);
  }

  std::vector<std::size_t> out_chains(out_offsets.back()),
      in_chains(in_offsets.back());
  auto out_c = out_offsets, in_c = in_offsets;

  // 5. Trace paths from end nodes, fill CSR directly.
  std::vector<bool> edge_used(sorted_edges.size());
  std::vector<std::size_t> idx, off{0};
  std::vector<int> powers;

  for (std::size_t i = 0; i < node_num; i++) {
    if (!is_end[i])
      continue;
    for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++) {
      if (edge_used[j])
        continue;
      std::size_t chain_idx = off.size() - 1;
      edge_used[j] = true;
      idx.push_back(i);
      powers.push_back(sorted_edges[j].power);
      std::size_t cur = sorted_edges[j].end;
      while (!is_end[cur]) {
        idx.push_back(cur);
        std::size_t nj = edge_offsets[cur];
        edge_used[nj] = true;
        cur = sorted_edges[nj].end;
      }
      idx.push_back(cur);
      off.push_back(idx.size());
      out_chains[out_c[compact_of[i]]++] = chain_idx;
      in_chains[in_c[compact_of[cur]]++] = chain_idx;
    }
  }

  // 6. Trace cycles from remaining edges, append directly to CSR.
  for (std::size_t i = 0; i < node_num; i++) {
    for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++) {
      if (edge_used[j])
        continue;
      std::size_t chain_idx = off.size() - 1;
      std::size_t start = i;
      edge_used[j] = true;
      idx.push_back(start);
      powers.push_back(sorted_edges[j].power);
      std::size_t cur = sorted_edges[j].end;
      while (cur != start) {
        idx.push_back(cur);
        std::size_t nj = edge_offsets[cur];
        edge_used[nj] = true;
        cur = sorted_edges[nj].end;
      }
      idx.push_back(cur);
      off.push_back(idx.size());

      out_offsets.push_back(out_offsets.back() + 1);
      in_offsets.push_back(in_offsets.back() + 1);
      out_chains.push_back(chain_idx);
      in_chains.push_back(chain_idx);
    }
  }

  return {std::move(idx),         std::move(off),        std::move(powers),
          std::move(out_offsets), std::move(out_chains), std::move(in_offsets),
          std::move(in_chains)};
}

} // namespace best_clipper
