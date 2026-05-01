#include "basic.hpp"
#include <iostream>

int main() {
    // Simple case: two overlapping squares
    multi_polygon first, second;
    bg::read_wkt("MULTIPOLYGON(((0 0, 0 2, 2 2, 2 0, 0 0)))", first);
    bg::read_wkt("MULTIPOLYGON(((1 1, 1 3, 3 3, 3 1, 1 1)))", second);

    // Collect segments
    std::vector<segment> segs;
    bg::for_each_segment(first, [&](const auto& seg) {
        segment s;
        bg::set<0,0>(s,bg::get<0,0>(seg)); bg::set<0,1>(s,bg::get<0,1>(seg));
        bg::set<1,0>(s,bg::get<1,0>(seg)); bg::set<1,1>(s,bg::get<1,1>(seg));
        segs.push_back(s);
    });
    bg::for_each_segment(second, [&](const auto& seg) {
        segment s;
        bg::set<0,0>(s,bg::get<0,0>(seg)); bg::set<0,1>(s,bg::get<0,1>(seg));
        bg::set<1,0>(s,bg::get<1,0>(seg)); bg::set<1,1>(s,bg::get<1,1>(seg));
        segs.push_back(s);
    });

    std::cout << "Num segments: " << segs.size() << std::endl;

    auto [hot_pixels, edges_with_power] = construct_edges_with_power(segs);
    std::cout << "Num hot pixels: " << hot_pixels.size() << std::endl;
    std::cout << "Num edges: " << edges_with_power.size() << std::endl;

    for (size_t i = 0; i < hot_pixels.size(); i++) {
        std::cout << "  hp[" << i << "] = " << bg::wkt(hot_pixels[i]) << std::endl;
    }
    for (auto& e : edges_with_power) {
        std::cout << "  edge: " << e.start << " -> " << e.end << "  power=" << e.power << std::endl;
    }

    // Build chains
    auto sorted_edges = edges_with_power;
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
    {
        std::vector<std::size_t> cnt(hot_pixels.size(), 0);
        for (auto& e : sorted_edges) {
            auto ts = e.power > 0 ? e.start : e.end;
            cnt[ts]++;
        }
        std::size_t cur = 0;
        for (std::size_t v = 0; v < hot_pixels.size(); v++) {
            edge_offsets[v] = cur;
            cur += cnt[v];
        }
        edge_offsets.back() = cur;
    }

    std::cout << "Edge offsets: ";
    for (auto o : edge_offsets) std::cout << o << " ";
    std::cout << std::endl;

    auto chains = build_chains(sorted_edges, edge_offsets, hot_pixels.size());
    std::cout << "Num chains: " << (chains.offsets.size() - 1) << std::endl;
    for (size_t c = 0; c + 1 < chains.offsets.size(); c++) {
        std::cout << "  chain " << c << " [p=" << chains.powers[c] << "]: ";
        for (size_t k = chains.offsets[c]; k < chains.offsets[c+1]; k++) {
            std::cout << chains.indices[k];
            if (k+1 < chains.offsets[c+1]) std::cout << " -> ";
        }
        std::cout << std::endl;
    }

    // Build face graph
    auto fg = build_face_graph(chains, hot_pixels);
    std::cout << "Num dcs: " << fg.sorted_dcs.size() << std::endl;

    // Dump face graph
    for (size_t v = 0; v < hot_pixels.size(); v++) {
        auto beg = fg.dc_begin[v];
        auto end = fg.dc_end[v];
        if (beg == end) continue;
        std::cout << "  vertex " << v << " (" << bg::wkt(hot_pixels[v]) << "): ";
        for (auto i = beg; i < end; i++) {
            auto dc = fg.sorted_dcs[i];
            std::cout << "dc" << dc.id << "(c" << dc.chain_id()
                      << (dc.is_forward()?"f":"r")
                      << " src=" << dc.source_node(chains)
                      << " tgt=" << dc.target_node(chains)
                      << " next_src=" << dc.next_along_source(chains)
                      << " p=" << dc.power(chains) << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << "Coplanar pairs:" << std::endl;
    for (auto [a,b] : fg.coplanar_pairs) {
        std::cout << "  dc" << a << " ~ dc" << b << std::endl;
    }

    // Compute winding numbers
    auto winding = compute_winding_numbers(chains, hot_pixels, fg);
    std::cout << "Num faces: " << winding.num_faces << std::endl;
    std::cout << "DC face winding numbers:" << std::endl;
    for (size_t i = 0; i < winding.dc_face_winding.size(); i++) {
        std::cout << "  dc" << i << ": face_winding=" << winding.dc_face_winding[i] << std::endl;
    }

    auto filter_union = [](int w) { return w > 0; };
    auto filter_inter = [](int w) { return w > 1; };

    std::cout << "\nUnion filter:" << std::endl;
    for (size_t i = 0; i < winding.dc_face_winding.size(); i++) {
        bool survive = filter_union(winding.dc_face_winding[i]);
        std::cout << "  dc" << i << ": survive=" << survive << std::endl;
    }

    std::cout << "\nIntersection filter:" << std::endl;
    for (size_t i = 0; i < winding.dc_face_winding.size(); i++) {
        bool survive = filter_inter(winding.dc_face_winding[i]);
        std::cout << "  dc" << i << ": survive=" << survive << std::endl;
    }

    return 0;
}
