#include <cstdio>
#include "core.hpp"
using namespace best_clipper;

int main() {
    multi_polygon a, b;
    bg::read_wkt("MULTIPOLYGON(((1941 2867,1964 2492,1818 2486,1941 2867)))", a);
    bg::read_wkt("MULTIPOLYGON(((1780 2877,1946 2821,1598 2541,1192 2638,1780 2877)))", b);

    auto segs = collect_segments(a, b);

    // Replicate construct_graph and check seg[1] specifically
    std::vector<std::pair<box, std::size_t>> boxes(segs.size());
    for (std::size_t i = 0; i < segs.size(); i++)
        boxes[i] = {bg::return_envelope<box>(segs[i]), i};

    bg::index::rtree<std::pair<box, std::size_t>, bg::index::quadratic<128>> rtree(std::move(boxes));

    std::vector<point> hot_pixels;
    for (const auto& segment : segs) {
        hot_pixels.emplace_back(bg::get<0, 0>(segment), bg::get<0, 1>(segment));
        hot_pixels.emplace_back(bg::get<1, 0>(segment), bg::get<1, 1>(segment));
    }

    for (std::size_t i = 0; i < segs.size(); i++)
        std::for_each(
            rtree.qbegin(bg::index::intersects(boxes[i].first)), rtree.qend(),
            [&](auto const& o) {
                if (auto p = get_intersection(segs[i], segs[o.second]))
                    hot_pixels.push_back(p.value());
            });

    std::sort(std::begin(hot_pixels), std::end(hot_pixels), [](auto p1, auto p2) {
        return std::pair{bg::get<0>(p1), bg::get<1>(p1)} <
               std::pair{bg::get<0>(p2), bg::get<1>(p2)};
    });
    auto last = std::unique(std::begin(hot_pixels), std::end(hot_pixels),
                            [](auto p1, auto p2) { return bg::equals(p1, p2); });
    hot_pixels.erase(last, hot_pixels.end());

    fprintf(stderr, "hot_pixels after dedup: %zu\n", hot_pixels.size());
    for (size_t i = 0; i < hot_pixels.size(); i++)
        fprintf(stderr, "  hp[%zu]: (%u,%u)\n", i, bg::get<0>(hot_pixels[i]), bg::get<1>(hot_pixels[i]));

    // Check: for segment 1 ((1964,2492)-(1818,2486)), which hot pixels lie on it?
    segment seg1 = segs[1];
    fprintf(stderr, "\nSegment 1: (%u,%u)-(%u,%u)\n",
        bg::get<0,0>(seg1), bg::get<0,1>(seg1),
        bg::get<1,0>(seg1), bg::get<1,1>(seg1));

    for (size_t i = 0; i < hot_pixels.size(); i++) {
        bool on = is_point_on_segment(hot_pixels[i], seg1);
        if (on)
            fprintf(stderr, "  hp[%zu]=(%u,%u) ON seg1\n", i, bg::get<0>(hot_pixels[i]), bg::get<1>(hot_pixels[i]));
    }

    // Check: rtree query for expanded pixels
    fprintf(stderr, "\nrtree query for seg[1] endpoints:\n");
    for (size_t i = 0; i < hot_pixels.size(); i++) {
        uint32_t x = bg::get<0>(hot_pixels[i]), y = bg::get<1>(hot_pixels[i]);
        // Only check points that could be seg1 endpoints
        if (bg::equals(hot_pixels[i], point{1964,2492}) || bg::equals(hot_pixels[i], point{1818,2486})) {
            auto min_corner = point{x > 0 ? x - 1 : 0, y > 0 ? y - 1 : 0};
            auto max_corner = point{x < UINT32_MAX ? x + 1 : UINT32_MAX, y < UINT32_MAX ? y + 1 : UINT32_MAX};
            fprintf(stderr, "  hp[%zu]=(%u,%u) expanded box: ((%u,%u),(%u,%u))\n",
                i, x, y,
                bg::get<0>(min_corner), bg::get<1>(min_corner),
                bg::get<0>(max_corner), bg::get<1>(max_corner));
            int count = 0;
            std::for_each(
                rtree.qbegin(bg::index::intersects(box{min_corner, max_corner})), rtree.qend(),
                [&](auto const& val) {
                    count++;
                    fprintf(stderr, "    rtree found: seg[%zu] box ((%u,%u),(%u,%u))\n",
                        val.second,
                        bg::get<0,0>(val.first), bg::get<0,1>(val.first),
                        bg::get<1,0>(val.first), bg::get<1,1>(val.first));
                });
            if (count == 0) fprintf(stderr, "    NO SEGMENTS FOUND!\n");
        }
    }

    return 0;
}
