#include "graph_helper.hpp"

namespace best_clipper {

std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>>& edges) {

    std::vector<std::size_t> parent(n);
    for (std::size_t i = 0; i < n; i++) parent[i] = i;

    auto find = [&](std::size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    for (auto [a, b] : edges) {
        auto ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    }

    std::vector<std::size_t> comp(n);
    std::vector<std::size_t> root_to_id(n, ~0ULL);
    std::size_t num_components = 0;
    for (std::size_t v = 0; v < n; v++) {
        auto r = find(v);
        if (root_to_id[r] == ~0ULL) root_to_id[r] = num_components++;
        comp[v] = root_to_id[r];
    }

    return comp;
}

}