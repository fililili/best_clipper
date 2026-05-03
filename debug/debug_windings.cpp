#include "basic.hpp"
#include <iostream>

int main() {
    multi_polygon first, second;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1), (6 6, 8 6, 8 8, 6 8, 6 6)))", first);
    bg::read_wkt("MULTIPOLYGON(((2 2, 2 7, 7 7, 7 2, 2 2)))", second);

    auto segs = collect_segments(first, second);
    auto [hot_pixels, ewp] = construct_edges_with_power(segs);

    std::cerr << "Hot pixels: " << hot_pixels.size() << "\n";
    for (size_t i = 0; i < hot_pixels.size(); i++)
        std::cerr << "  hp[" << i << "]=" << bg::wkt(hot_pixels[i]) << "\n";

    std::cerr << "Edges with power: " << ewp.size() << "\n";
    for (auto& e : ewp)
        std::cerr << "  " << e.start << "->" << e.end << " p=" << e.power << "\n";

    // Sort edges for chain building
    auto sorted_edges = ewp;
    std::sort(sorted_edges.begin(), sorted_edges.end(),
              [](const edge_with_power_t& a, const edge_with_power_t& b) {
                  auto ts = [](const edge_with_power_t& e) { return e.power > 0 ? e.start : e.end; };
                  auto sa = ts(a), sb = ts(b);
                  if (sa != sb) return sa < sb;
                  auto ta = a.power > 0 ? a.end : a.start;
                  auto tb = b.power > 0 ? b.end : b.start;
                  return ta < tb;
              });
    std::vector<std::size_t> edge_offsets(hot_pixels.size() + 1, 0);
    { std::vector<std::size_t> cnt(hot_pixels.size(), 0);
      for (auto& e : sorted_edges) cnt[e.power > 0 ? e.start : e.end]++;
      std::size_t cur = 0;
      for (std::size_t v = 0; v < hot_pixels.size(); v++) { edge_offsets[v] = cur; cur += cnt[v]; }
      edge_offsets.back() = cur; }

    auto chains = build_chains(sorted_edges, edge_offsets, hot_pixels.size());
    std::cerr << "Chains: " << (chains.offsets.size()-1) << "\n";
    for (size_t c = 0; c + 1 < chains.offsets.size(); c++) {
        std::cerr << "  chain[" << c << "] p=" << chains.powers[c] << ": ";
        for (size_t k = chains.offsets[c]; k < chains.offsets[c+1]; k++)
            std::cerr << chains.indices[k] << " ";
        std::cerr << "\n";
    }

    auto [sorted_hcs, hc_begin, hc_end, next_hc, prev_hc, coplanar] = build_hc_graph(chains, hot_pixels);
    std::cerr << "HCs: " << (chains.offsets.size()-1)*2 << "\n";
    for (size_t v = 0; v < hot_pixels.size(); v++) {
        auto beg = hc_begin[v], end = hc_end[v];
        if (beg == end) continue;
        std::cerr << "  vertex " << v << " (" << bg::wkt(hot_pixels[v]) << "): ";
        for (auto i = beg; i < end; i++) {
            auto hc = sorted_hcs[i];
            std::cerr << "hc" << hc.id << "(" << (hc.is_forward()?"f":"r") << "c" << hc.chain_id()
                      << " src=" << hc.source_node(chains) << " nxt=" << hc.next_along_source(chains)
                      << " p=" << hc.power(chains) << ") ";
        }
        std::cerr << "\n";
    }

    std::cerr << "Coplanar pairs:\n";
    for (auto [a,b] : coplanar)
        std::cerr << "  hc" << a << " ~ hc" << b << "\n";

    std::cerr << "Next:\n";
    for (size_t i = 0; i < (chains.offsets.size()-1)*2; i++)
        std::cerr << "  hc" << i << " -> hc" << next_hc[i].id << "\n";

    // Debug vertex components
    {
        std::size_t nv = hot_pixels.size();
        std::vector<std::size_t> vp(nv);
        for (std::size_t i = 0; i < nv; i++) vp[i] = i;
        auto vfind2 = [&](std::size_t x) {
            std::size_t r = x;
            while (vp[r] != r) r = vp[r];
            while (x != r) { std::size_t nxt = vp[x]; vp[x] = r; x = nxt; }
            return r;
        };
        for (auto hc : sorted_hcs) {
            std::size_t a = vfind2(hc.source_node(chains));
            std::size_t b = vfind2(hc.target_node(chains));
            if (a != b) vp[a] = b;
        }
        for (std::size_t v = 0; v < nv; v++) vp[v] = vfind2(v);
        std::cerr << "Vertex groups: ";
        for (std::size_t v = 0; v < nv; v++)
            std::cerr << "v" << v << "->g" << vp[v] << " ";
        std::cerr << "\n";
    }

    auto [ext_hcs, ray_pairs] = find_exterior(chains, hot_pixels, sorted_hcs, hc_begin, hc_end);
    std::cerr << "Exterior HCs: ";
    for (auto e : ext_hcs) std::cerr << "hc" << e << " ";
    std::cerr << "\n";
    std::cerr << "Ray coplanar pairs: ";
    for (auto [a,b] : ray_pairs) std::cerr << "(" << a << "," << b << ") ";
    std::cerr << "\n";

    auto [hc_winding, face_root] = compute_hc_winding(chains, sorted_hcs, hc_begin, hc_end, coplanar, ray_pairs, ext_hcs);
    std::cerr << "HC windings:\n";
    for (size_t i = 0; i < hc_winding.size(); i++)
        std::cerr << "  hc" << i << " w=" << hc_winding[i] << " face_root=" << face_root[i] << "\n";

    std::cerr << "Survive union (w>0): ";
    for (size_t i = 0; i < hc_winding.size(); i++)
        if (hc_winding[i] > 0) std::cerr << "hc" << i << " ";
    std::cerr << "\n";
    std::cerr << "Survive intersection (w>1): ";
    for (size_t i = 0; i < hc_winding.size(); i++)
        if (hc_winding[i] > 1) std::cerr << "hc" << i << " ";
    std::cerr << "\n";

    // Show face groups
    std::cerr << "Face groups:\n";
    for (size_t i = 0; i < face_root.size(); i++) {
        if (face_root[i] == i) {
            std::cerr << "  Face " << i << " (w=" << hc_winding[i] << "): ";
            for (size_t j = 0; j < face_root.size(); j++)
                if (face_root[j] == i) std::cerr << "hc" << j << " ";
            std::cerr << "\n";
        }
    }

    auto result_inter = run_pipeline(segs, [](int w) { return w > 1; });
    std::cerr << "Intersection area: " << bg::area(result_inter) << "\n";
    std::cerr << "Intersection num polygons: " << result_inter.size() << "\n";
    std::cerr << "Intersection valid: " << bg::is_valid(result_inter) << "\n";
    std::cerr << "Intersection WKT: " << bg::wkt(result_inter) << "\n";

    auto result_union = run_pipeline(segs, [](int w) { return w > 0; });
    std::cerr << "Union area: " << bg::area(result_union) << "\n";
    std::cerr << "Union num polygons: " << result_union.size() << "\n";
    std::cerr << "Union valid: " << bg::is_valid(result_union) << "\n";
    std::cerr << "Union WKT: " << bg::wkt(result_union) << "\n";
    if (result_union.size() > 0) {
        auto& p = result_union[0];
        std::cerr << "  outer area: " << bg::area(p.outer()) << "\n";
        std::cerr << "  num holes: " << p.inners().size() << "\n";
        for (size_t i = 0; i < p.inners().size(); i++) {
            std::cerr << "  hole[" << i << "] area: " << bg::area(p.inners()[i]) << "\n";
        }
        std::cerr << "  polygon valid: " << bg::is_valid(p) << "\n";
    }

    // Check expected union
    multi_polygon expected_union;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 2, 2 2, 2 3, 1 3, 1 1), (8 8, 6 8, 6 7, 7 7, 7 6, 8 6, 8 8)))", expected_union);
    std::cerr << "Expected union area: " << bg::area(expected_union) << std::endl;
    std::cerr << "Expected union valid: " << bg::is_valid(expected_union) << std::endl;
    std::cerr << "Union valid: " << bg::is_valid(result_union) << std::endl;
    std::cerr << "Equals expected union: " << bg::equals(result_union, expected_union) << "\n";
    return 0;
}
