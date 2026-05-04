#include "core.hpp"
#include <iostream>
using namespace best_clipper;

int main() {
    // Simple case: two overlapping squares
    multi_polygon first, second;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 2, 2 2, 2 0, 0 0)))", first);
    bg::read_wkt("MULTIPOLYGON(((1 1, 1 3, 3 3, 3 1, 1 1)))", second);

    auto segs = collect_segments(first, second);
    std::cout << "Num segments: " << segs.size() << std::endl;

    auto [hot_pixels, edges_with_power] = construct_edges_with_power(segs);
    std::cout << "Num hot pixels: " << hot_pixels.size() << std::endl;
    std::cout << "Num edges: " << edges_with_power.size() << std::endl;

    for (size_t i = 0; i < hot_pixels.size(); i++)
        std::cout << "  hp[" << i << "] = " << bg::wkt(hot_pixels[i]) << std::endl;
    for (auto& e : edges_with_power)
        std::cout << "  edge: " << e.start << " -> " << e.end << "  power=" << e.power << std::endl;

    // Build chains (current API: 2 params, does internal sorting)
    auto chains = build_chains(edges_with_power, hot_pixels.size());
    std::cout << "Num chains: " << (chains.offsets.size() - 1) << std::endl;
    for (size_t c = 0; c + 1 < chains.offsets.size(); c++) {
        std::cout << "  chain " << c << " [p=" << chains.powers[c] << "]: ";
        for (size_t k = chains.offsets[c]; k < chains.offsets[c + 1]; k++) {
            std::cout << chains.indices[k];
            if (k + 1 < chains.offsets[c + 1]) std::cout << " -> ";
        }
        std::cout << std::endl;
    }

    // Build half-chain graph (current API)
    auto [sorted_hcs, hc_begin, hc_end, next_hc, coplanar] = build_half_chain_graph(chains, hot_pixels);
    std::cout << "Num HCs: " << sorted_hcs.size() << std::endl;

    // Dump half-chain graph
    for (size_t v = 0; v < hot_pixels.size(); v++) {
        auto beg = hc_begin[v], end = hc_end[v];
        if (beg == end) continue;
        std::cout << "  vertex " << v << " (" << bg::wkt(hot_pixels[v]) << "): ";
        for (auto i = beg; i < end; i++) {
            auto hc = sorted_hcs[i];
            std::cout << "hc" << hc.id << "(c" << hc.chain_id()
                      << (hc.is_forward() ? "f" : "r")
                      << " src=" << hc.source_node(chains)
                      << " tgt=" << hc.target_node(chains)
                      << " next_src=" << hc.next_along_source(chains)
                      << " p=" << hc.power(chains) << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << "Coplanar pairs:" << std::endl;
    for (auto [a, b] : coplanar)
        std::cout << "  hc" << a << " ~ hc" << b << std::endl;

    // Ray casting + winding
    auto [exterior_hcs, ray_pairs] =
        find_exterior(chains, hot_pixels, sorted_hcs, hc_begin, hc_end);
    std::cout << "Exterior HCs: ";
    for (auto e : exterior_hcs) std::cout << "hc" << e << " ";
    std::cout << std::endl;

    auto winding = compute_winding(chains, coplanar, ray_pairs, exterior_hcs);
    std::cout << "HC windings:" << std::endl;
    for (size_t i = 0; i < winding.size(); i++)
        std::cout << "  hc" << i << ": w=" << winding[i] << std::endl;

    auto filter_union = [](int w) { return w > 0; };
    auto filter_inter = [](int w) { return w > 1; };

    std::cout << "\nUnion filter:" << std::endl;
    for (size_t i = 0; i < winding.size(); i++)
        std::cout << "  hc" << i << ": survive=" << filter_union(winding[i]) << std::endl;

    std::cout << "\nIntersection filter:" << std::endl;
    for (size_t i = 0; i < winding.size(); i++)
        std::cout << "  hc" << i << ": survive=" << filter_inter(winding[i]) << std::endl;

    return 0;
}
