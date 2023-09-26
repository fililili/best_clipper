#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <ranges>
#include <optional>
#include <stack>
#include <random>
#include <chrono>
#include <string>
#include <numeric>
#include <boost/geometry.hpp>
#include <boost/container/flat_map.hpp>

//#define LOG_TIME

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<int>;
using segment = bg::model::segment<point>;
using box = bg::model::box<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;

enum class fake_bool {
    fake_false,
    fake_true
};

std::string get_time_string(auto t) {
    auto ti = std::chrono::system_clock::to_time_t(t);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&ti), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

auto log_done_time(std::string work) {
#ifdef LOG_TIME
    using namespace std::chrono_literals;
    static auto before_time = std::chrono::system_clock::now();
    auto after_time = std::chrono::system_clock::now();
    std::cout << work << " is done, before time: " << get_time_string(before_time) << ", after time: " << get_time_string(after_time) << ", runtime: " << (after_time - before_time) / 1s << "s" << std::endl;
    before_time = after_time;
#endif
}


auto construct_multi_polygons(auto&& rings) {
    log_done_time("start of construct multi polygons");
    std::vector<fake_bool> is_rings_cw(rings.size());
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (bg::area(rings[i]) > 0)
            is_rings_cw[i] = fake_bool::fake_true;
        else
            is_rings_cw[i] = fake_bool::fake_false;
    }
    std::vector<std::pair<box, std::size_t>> cw_rings_box;
    multi_polygon ret;
    log_done_time("judge cw");

    cw_rings_box.reserve(rings.size());
    ret.reserve(rings.size());
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (is_rings_cw[i] == fake_bool::fake_true) {
            cw_rings_box.emplace_back(bg::return_envelope<box>(rings[i]), ret.size());
            polygon poly{ std::move(rings[i]) };
            ret.emplace_back(std::move(poly));
        }
    }
    log_done_time("build multi polygon");

    bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > rings_box_tree(cw_rings_box);
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (is_rings_cw[i] == fake_bool::fake_false) {
            for (auto itr = rings_box_tree.qbegin(bg::index::contains(rings[i][0])); itr != rings_box_tree.qend(); ++itr) {
                if (bg::relate(rings[i][0], bg::exterior_ring(ret[itr->second]), bg::de9im::static_mask<'T', 'F', 'F'>{})) {
                    bg::interior_rings(ret[itr->second]).emplace_back(std::move(rings[i]));
                    break;
                }
            }
        }
    }
    log_done_time("find hole parent");
    return ret;

}


// may be overflow, need study more
// current algorithm from https://leetcode.cn/problems/intersection-lcci/solutions/197813/jiao-dian-by-leetcode-solution/
std::optional<point> get_intersection(segment s1, segment s2) {
    using ctype = decltype(bg::get<0, 0>(s1));
    double x1 = bg::get<0, 0>(s1);
    double y1 = bg::get<0, 1>(s1);
    double x2 = bg::get<1, 0>(s1);
    double y2 = bg::get<1, 1>(s1);
    double x3 = bg::get<0, 0>(s2);
    double y3 = bg::get<0, 1>(s2);
    double x4 = bg::get<1, 0>(s2);
    double y4 = bg::get<1, 1>(s2);
    if ((y4 - y3) * (x2 - x1) == (y2 - y1) * (x4 - x3)) return {}; // return empty when paraller
    double t1 = (double)(x3 * (y4 - y3) + y1 * (x4 - x3) - y3 * (x4 - x3) - x1 * (y4 - y3)) / ((x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1));
    double t2 = (double)(x1 * (y2 - y1) + y3 * (x2 - x1) - y1 * (x2 - x1) - x3 * (y2 - y1)) / ((x4 - x3) * (y2 - y1) - (x2 - x1) * (y4 - y3));
    // 判断 t1 和 t2 是否均在 (0, 1) 之间
    if (t1 > 0.0 && t1 < 1.0 && t2 > 0.0 && t2 < 1.0) {
        // -0.5 < x <= 0.5 to 0
        return point{ static_cast<ctype>(std::ceil(x1 + t1 * (x2 - x1) - 0.5)), static_cast<ctype>(std::ceil(y1 + t1 * (y2 - y1) - 0.5)) };
    }
    else {
        return {};
    }
}

bool is_point_on_segment(point p, segment s) {
    // will use better algorithm
    auto x = bg::get<0>(p);
    auto y = bg::get<1>(p);
    auto x1 = bg::get<0, 0>(s);
    auto y1 = bg::get<0, 1>(s);
    auto x2 = bg::get<1, 0>(s);
    auto y2 = bg::get<1, 1>(s);
    if ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) < 0 || (x2 - x1) * (x2 - x) + (y2 - y1) * (y2 - y) < 0) return false;

    if ((2 * x + 1 - x1 - x2) * (y1 - y2) > (2 * y + 1 - y1 - y2) * (x1 - x2) &&
        (2 * x + 1 - x1 - x2) * (y1 - y2) >= (2 * y - 1 - y1 - y2) * (x1 - x2) &&
        (2 * x - 1 - x1 - x2) * (y1 - y2) >= (2 * y + 1 - y1 - y2) * (x1 - x2) &&
        (2 * x - 1 - x1 - x2) * (y1 - y2) >= (2 * y - 1 - y1 - y2) * (x1 - x2)) {
        return false;
    }
    if ((2 * x + 1 - x1 - x2) * (y1 - y2) < (2 * y + 1 - y1 - y2) * (x1 - x2) &&
        (2 * x + 1 - x1 - x2) * (y1 - y2) <= (2 * y - 1 - y1 - y2) * (x1 - x2) &&
        (2 * x - 1 - x1 - x2) * (y1 - y2) <= (2 * y + 1 - y1 - y2) * (x1 - x2) &&
        (2 * x - 1 - x1 - x2) * (y1 - y2) <= (2 * y - 1 - y1 - y2) * (x1 - x2)) {
        return false;
    }
    return true;
}

// construct a graph with edge property (power) by segs
// segs should can create rings
auto construct_rings(const auto& segs, auto filter) {
    log_done_time("start of construct_rings");
    std::vector<std::pair<box, std::size_t>> boxes(segs.size());
    for (std::size_t i = 0; i < segs.size(); i++) {
        boxes[i] = { bg::return_envelope<box>(segs[i]), i };
    }
    bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > segs_box_rtree(boxes);

    std::vector<point> hot_pixels;
    hot_pixels.reserve(segs.size() * 2);

    for (const auto& seg : segs) {
        hot_pixels.emplace_back(bg::get<0, 0>(seg), bg::get<0, 1>(seg));
    }

    log_done_time("build segs rtree");
    for (std::size_t i = 0; i < segs.size(); i++) {
        std::for_each(segs_box_rtree.qbegin(bg::index::intersects(boxes[i].first)), segs_box_rtree.qend(),
            [&](auto const& other_seg) {
                std::optional<point> p = get_intersection(segs[i], segs[other_seg.second]); // find other_seg seg intersection
                if (p) {
                    hot_pixels.push_back(p.value());
                }
            }
        );

    }
    log_done_time("find hot_pixels");
    
    {
        std::sort(std::begin(hot_pixels), std::end(hot_pixels),
            [](auto p1, auto p2) {
                // can dig more to see if we can use a better sort algorithm to better order pixels
                return std::pair{ bg::get<0>(p1), bg::get<1>(p1) } < std::pair{ bg::get<0>(p2), bg::get<1>(p2) };
            }
        ); // need define a better one
        auto last = std::unique(std::begin(hot_pixels), std::end(hot_pixels), [](auto p1, auto p2) { return bg::equals(p1, p2); });
        hot_pixels.erase(last, hot_pixels.end());

        // reorder hot pixels to make memory cache better
        bg::index::rtree<point, bg::index::quadratic<128> > hot_pixels_rtree{ hot_pixels };
        hot_pixels = std::vector<point>{ std::begin(hot_pixels_rtree), std::end(hot_pixels_rtree) };
        
    }
    log_done_time("order hot pixels");

    std::vector<std::pair<std::size_t, std::size_t> > seg_pixel_pairs;
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        constexpr auto expand = 1;
        auto min_corner = point{ bg::get<0>(hot_pixels[i]) - expand, bg::get<1>(hot_pixels[i]) - expand };
        auto max_corner = point{ bg::get<0>(hot_pixels[i]) + expand, bg::get<1>(hot_pixels[i]) + expand };

        std::for_each(segs_box_rtree.qbegin(bg::index::intersects(box{ min_corner, max_corner })), segs_box_rtree.qend(),
            [&](auto const& val) {
                if (is_point_on_segment(hot_pixels[i], segs[val.second])) {
                    seg_pixel_pairs.emplace_back(val.second, i);
                }
            }
        );

    }
    log_done_time("find segs on hot_pixels");

    // can be more precise
    struct less_by_segment {
        bool operator() (std::pair<std::size_t, std::size_t>  v1, std::pair<std::size_t, std::size_t>  v2) {
            auto p1 = _hot_pixels[v1.second];
            auto p2 = _hot_pixels[v2.second];
            double x1 = bg::get<0>(p1);
            double y1 = bg::get<1>(p1);
            double x2 = bg::get<0>(p2);
            double y2 = bg::get<1>(p2);
            double x3 = bg::get<0, 0>(s);
            double y3 = bg::get<0, 1>(s);
            double x4 = bg::get<1, 0>(s);
            double y4 = bg::get<1, 1>(s);

            return (x2 - x1) * (x4 - x3) + (y2 - y1) * (y4 - y3) > 0;
        }
        segment s;
        const std::vector<point>& _hot_pixels;
    };
    std::vector<std::pair<std::size_t, std::size_t> > edges;
    {
        std::vector<int> segs_times(segs.size());
        for (auto [seg, _] : seg_pixel_pairs) {
            segs_times[seg]++;
        }
        std::vector<int> begin_location(segs_times.size());
        std::exclusive_scan(std::begin(segs_times), std::end(segs_times), std::begin(begin_location), 0);
        decltype(seg_pixel_pairs) _seg_pixel_pairs{seg_pixel_pairs.size()};
        auto current_location{begin_location};
        for (std::size_t i = 0; i < seg_pixel_pairs.size(); i++) {
            _seg_pixel_pairs[current_location[seg_pixel_pairs[i].first]++] = seg_pixel_pairs[i];
        }
        seg_pixel_pairs = std::move(_seg_pixel_pairs);
        auto end_location{current_location};
        for (std::size_t i = 0; i < segs.size(); i++) {
            auto begin = std::begin(seg_pixel_pairs);
            auto cur_begin = begin + begin_location[i];
            auto cur_end = begin + end_location[i];
            auto cur_last = cur_end - 1;
            std::sort(cur_begin, cur_end, less_by_segment{ segs[i], hot_pixels } );
            
            for (; cur_begin != cur_last; cur_begin++) {
                edges.emplace_back(cur_begin->second, std::next(cur_begin)->second);
            }
        }
    }
    log_done_time("build edges");

    constexpr auto to_order = [](auto e) {
        return decltype(e){ (std::max)(std::get<0>(e), std::get<1>(e)), (std::min)(std::get<0>(e), std::get<1>(e)) };
        };

    std::sort(std::begin(edges), std::end(edges), [&](auto e1, auto e2) {
        return to_order(e1) < to_order(e2);
        }
    );

    std::vector<int> edges_power(edges.size());
    // remove duplicated edges and calculate edges power
    {
        std::size_t i = 0;
        std::size_t j = 0;
        auto cur_first = 0;
        auto cur_power = 0;
        for (std::size_t i = 0; i < edges.size(); i++) {
            if (to_order(edges[i]) == edges[i]) {
                cur_power += 1;
            }
            else {
                cur_power -= 1;
            }
            if (i == edges.size() - 1 || to_order(edges[i + 1]) != to_order(edges[cur_first])) {
                if (cur_power == 0) {
                    cur_first = i + 1;
                }
                else {
                    edges[j] = to_order(edges[i]);
                    edges_power[j] = cur_power;
                    cur_power = 0;
                    j++;
                    cur_first = i + 1;
                }
            }
        }
        edges.erase(std::begin(edges) + j, std::end(edges));
        edges_power.erase(std::begin(edges_power) + j, std::end(edges_power));
    }

    log_done_time("build edges power");

    auto source = [&](auto de) {
        if (de.first) return edges[de.second].first;
        else return edges[de.second].second;
        };
    auto target = [&](auto de) {
        if (de.first) return edges[de.second].second;
        else return edges[de.second].first;
        };
    auto power = [&](auto de) {
        if (de.first) return edges_power[de.second];
        else return -edges_power[de.second];
        };

    auto get_direction = [&](auto i) {
        // can be more precise
        double dx = bg::get<0>(hot_pixels[target(i)]) - bg::get<0>(hot_pixels[source(i)]);
        double dy = bg::get<1>(hot_pixels[target(i)]) - bg::get<1>(hot_pixels[source(i)]);
        enum class quadrant { _1, _2, _3, _4, zero };
        if (dx > 0 && dy >= 0) {
            return std::pair{ quadrant::_1, dy / dx };
        }
        else if (dx <= 0 && dy > 0) {
            return std::pair{ quadrant::_2, -dx / dy };
        }
        else if (dx < 0 && dy <= 0) {
            return std::pair{ quadrant::_3, dy / dx };
        }
        else if (dx >= 0 && dy < 0) {
            return std::pair{ quadrant::_4, -dx / dy };
        }
        // never happen
        return std::pair{ quadrant::zero, 0.0 };
        };
        
    std::vector<std::pair<bool, std::size_t> > direct_edges(edges.size() * 2);
    {
        std::vector<int> hot_pixels_times(hot_pixels.size());
        for (auto edge : edges) {
            hot_pixels_times[edge.first]++;
            hot_pixels_times[edge.second]++;
        }
        std::vector<int> begin_location(hot_pixels.size());
        std::exclusive_scan(std::begin(hot_pixels_times), std::end(hot_pixels_times), std::begin(begin_location), 0);
        assert(edges.size() * 2 == hot_pixels_times.back() + begin_location.back());
        auto current_location{begin_location};
        for (std::size_t i = 0; i < edges.size(); i++) {
            direct_edges[current_location[edges[i].first]++] = { true, i };
            direct_edges[current_location[edges[i].second]++] = { false, i };
        }
        auto end_location{current_location};
        assert(edges.size() * 2 == end_location.back());
        for (std::size_t i = 0; i < hot_pixels.size(); i++) {
            if (end_location[i] - begin_location[i] < 3) continue;
            auto begin = std::begin(direct_edges);
            std::sort(begin + begin_location[i], begin + end_location[i],
                [&](auto i1, auto i2) {
                    return get_direction(i1) < get_direction(i2);
                }
            );
        }
    }
    log_done_time("sort direct edges");

    std::vector<std::pair<std::size_t, std::size_t> > direct_edge_pairs(edges.size());
    for (std::size_t i = 0; i < direct_edges.size(); i++) {
        if (direct_edges[i].first) direct_edge_pairs[direct_edges[i].second].second = i;
        else direct_edge_pairs[direct_edges[i].second].first = i;
    }
    auto get_dual = [&](auto i) {
        if (direct_edges[i].first) {
            return direct_edge_pairs[direct_edges[i].second].first;
        }
        else {
            return direct_edge_pairs[direct_edges[i].second].second;
        }
        };

    for (std::size_t i = 0; i < direct_edges.size(); i++) {
        auto dual = get_dual(i);
        //std::cout << "i = " << i << ", dual i = " << dual << std::endl;
    }

    // link direct edges, except the outer face, all faces are cw oriented
    std::vector<std::size_t> next_direct_edges(direct_edges.size());
    std::vector<std::size_t> pre_direct_edges(direct_edges.size());
    {
        std::size_t cur_first = 0;
        for (std::size_t i = 0; i < direct_edges.size(); i++) {
            //std::cout << "source: " << source(direct_edges[i]) << ":  " << bg::wkt(hot_pixels[source(direct_edges[i])]) << "  ";
            //std::cout << "target: " << target(direct_edges[i]) << ":  " << bg::wkt(hot_pixels[target(direct_edges[i])]) << "  ";
            //std::cout << "power: " << power(direct_edges[i]) << std::endl;
            if (i + 1 == direct_edges.size() || source(direct_edges[i + 1]) != source(direct_edges[cur_first])) {
                pre_direct_edges[cur_first] = get_dual(i);
                next_direct_edges[get_dual(i)] = cur_first;
                cur_first = i + 1;
            }
            else {
                pre_direct_edges[i + 1] = get_dual(i);
                next_direct_edges[get_dual(i)] = i + 1;
            }
        }
    }
    log_done_time("connect direct edges");

    std::vector<std::size_t> edges_face_id(direct_edges.size());
    std::vector<std::pair<ring, std::size_t> > cw_faces;
    std::size_t cur_face_id = 1;
    boost::container::flat_map<std::size_t, std::size_t> ccw_face_id_to_direct_edge_id;
    for (std::size_t i = 0; i < direct_edges.size(); i++) {
        if (edges_face_id[i] != 0) continue;
        auto cur_first_id = i;
        ring r;
        // ring is closed, need to push one duplicated point
        r.push_back(hot_pixels[target(direct_edges[cur_first_id])]);
        do {
            edges_face_id[i] = cur_face_id;
            i = next_direct_edges[i];
            r.push_back(hot_pixels[target(direct_edges[i])]);
        } while (cur_first_id != i);
        if (bg::area(r) > 0) { // for cw orientation, the area > 0, may use more precise algorithm
            cw_faces.emplace_back(r, cur_face_id);
        }
        else {
            ccw_face_id_to_direct_edge_id.insert(std::pair{ cur_face_id, cur_first_id });
        }
        cur_face_id++;
    }

    log_done_time("build faces");

    boost::container::flat_multimap<std::size_t, std::size_t> face_contain_relations;
    {
        std::vector<std::pair<box, std::size_t> > cw_faces_box(cw_faces.size());
        for (std::size_t i = 0; i < cw_faces.size(); i++) {
            cw_faces_box[i].first = bg::return_envelope<box>(cw_faces[i].first);
            cw_faces_box[i].second = i;
        }
        bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > faces_rtree{ cw_faces_box };
        for (auto [face_id, direct_edge_id] : ccw_face_id_to_direct_edge_id) {
            point p = hot_pixels[source(direct_edges[direct_edge_id])];
            auto itr = faces_rtree.qbegin(bg::index::contains(p));
            if (itr == faces_rtree.qend())
                face_contain_relations.insert(std::pair<std::size_t, std::size_t>{ 0, face_id });
            else {
                for (; itr != faces_rtree.qend(); ++itr) {
                    if (bg::relate(p, cw_faces[itr->second].first, bg::de9im::static_mask<'T', 'F', 'F'>{})) {
                        face_contain_relations.insert(std::pair<std::size_t, std::size_t>{ cw_faces[itr->second].second , face_id });
                        break;
                    }
                }
            }
        }
    }
    log_done_time("build face contain relations");

    std::vector<std::optional<int> > faces_cw_power(cur_face_id);
    std::vector<fake_bool> direct_edges_exist(direct_edges.size());
    std::stack<std::pair<std::size_t, std::size_t> > stk;
    faces_cw_power[0] = 0;
    for (auto [b, e] = face_contain_relations.equal_range(0); b != e; b++) {
        auto next_face_id = (*b).second;
        if (!faces_cw_power[next_face_id]) {
            faces_cw_power[next_face_id] = 0;
            stk.push({ next_face_id, ccw_face_id_to_direct_edge_id[next_face_id] });
        }
    }
    while (!stk.empty()) {
        auto [cur_face_id, cur_first_direct_edge] = stk.top();
        stk.pop();
        for (auto [b, e] = face_contain_relations.equal_range(cur_face_id); b != e; b++) {
            auto next_face_id = (*b).second;
            if (!faces_cw_power[next_face_id]) {
                faces_cw_power[next_face_id] = faces_cw_power[cur_face_id];
                stk.push({ next_face_id, ccw_face_id_to_direct_edge_id[next_face_id] });
            }
        }

        auto cur_direct_edge = cur_first_direct_edge;
        do {
            if (filter(faces_cw_power[cur_face_id].value())) direct_edges_exist[cur_direct_edge] = fake_bool::fake_true;
            else direct_edges_exist[cur_direct_edge] = fake_bool::fake_false;

            auto dual = get_dual(cur_direct_edge);
            auto next_face_id = edges_face_id[dual];
            if (!faces_cw_power[next_face_id]) {
                faces_cw_power[next_face_id] = faces_cw_power[cur_face_id].value() + power(direct_edges[dual]);
                stk.push({ next_face_id, dual });
            }
            cur_direct_edge = next_direct_edges[cur_direct_edge];
        } while (cur_direct_edge != cur_first_direct_edge);
    }
    log_done_time("traversal faces");

    for (auto&& [de1, de2] : direct_edge_pairs) {
        if (direct_edges_exist[de1] == direct_edges_exist[de2]) {
            direct_edges_exist[de1] = fake_bool::fake_false;
            direct_edges_exist[de2] = fake_bool::fake_false;
            auto pre1 = pre_direct_edges[de1];
            auto pre2 = pre_direct_edges[de2];
            auto next1 = next_direct_edges[de1];
            auto next2 = next_direct_edges[de2];
            pre_direct_edges[next1] = pre2;
            pre_direct_edges[next2] = pre1;
            next_direct_edges[pre1] = next2;
            next_direct_edges[pre2] = next1;
            // no use in fact
            pre_direct_edges[de1] = de2;
            pre_direct_edges[de2] = de1;
            next_direct_edges[de1] = de2;
            next_direct_edges[de2] = de1;
        }
    }

    std::vector<ring> ret_rings;
    {
        
        for (std::size_t i = 0; i < direct_edges_exist.size(); i++) {

            if (direct_edges_exist[i] == fake_bool::fake_false) continue;

            std::size_t cur_first_id = i;
            ring r;
            // ring is closed, need to push one duplicated point
            r.push_back(hot_pixels[target(direct_edges[i])]);
            direct_edges_exist[i] = fake_bool::fake_false;
            do {
                i = next_direct_edges[i];
                r.push_back(hot_pixels[target(direct_edges[i])]);
                direct_edges_exist[i] = fake_bool::fake_false;
            } while (cur_first_id != i);
            ret_rings.push_back(std::move(r));
            cur_first_id = i;
        }
    }
    for (auto&& ring : ret_rings) {
        //std::cout << bg::wkt(ring) << std::endl;
        //std::cout << "is valid: " << bg::is_valid(ring) << std::endl;
    }

    return ret_rings; // graph
}

auto add (const auto& ps1, const auto& ps2) {
    std::vector<segment > segs;
    segs.reserve(bg::num_segments(ps1) + bg::num_segments(ps2));
    bg::for_each_segment(ps1,
        [&](const auto& seg) {
            segment s;
            boost::geometry::set<0, 0>(s, boost::geometry::get<0, 0>(seg));
            boost::geometry::set<0, 1>(s, boost::geometry::get<0, 1>(seg));
            boost::geometry::set<1, 0>(s, boost::geometry::get<1, 0>(seg));
            boost::geometry::set<1, 1>(s, boost::geometry::get<1, 1>(seg));
            segs.emplace_back(s);
        }
    );
    bg::for_each_segment(ps2,
        [&](const auto& seg) {
            segment s;
            boost::geometry::set<0, 0>(s, boost::geometry::get<0, 0>(seg));
            boost::geometry::set<0, 1>(s, boost::geometry::get<0, 1>(seg));
            boost::geometry::set<1, 0>(s, boost::geometry::get<1, 0>(seg));
            boost::geometry::set<1, 1>(s, boost::geometry::get<1, 1>(seg));
            segs.emplace_back(s);
        }
    );

    constexpr auto filter = [](auto cw_power) { return cw_power > 0; };
    return construct_multi_polygons(construct_rings(segs, filter));
}
auto self_or(auto r) {
    std::vector<segment > segs;
    bg::for_each_segment(r,
        [&](const auto& seg) {
            segment s;
            boost::geometry::set<0, 0>(s, boost::geometry::get<0, 0>(seg));
            boost::geometry::set<0, 1>(s, boost::geometry::get<0, 1>(seg));
            boost::geometry::set<1, 0>(s, boost::geometry::get<1, 0>(seg));
            boost::geometry::set<1, 1>(s, boost::geometry::get<1, 1>(seg));
            segs.emplace_back(s);
        }
    );
    auto rings = construct_rings(segs, [](auto cw_power) { return cw_power > 0; });
    auto ret = construct_multi_polygons(std::move(rings));
    return ret;
}

auto benchmark(int size) {
    using namespace std::chrono_literals;
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(-600, 600);

    ring r1, r2;
    point p{ 0, 0 };
    r1.emplace_back(p);
    for (int i = 0; i < size; i++) {
        bg::set<0>(p, bg::get<0>(p) + distrib(gen));
        bg::set<1>(p, bg::get<1>(p) + distrib(gen));
        r1.emplace_back(p);
    }
    r2.emplace_back(0, 0);
    p = point{ 0, 0 };
    r2.emplace_back(p);
    for (int i = 0; i < size; i++) {
        bg::set<0>(p, bg::get<0>(p) + distrib(gen));
        bg::set<1>(p, bg::get<1>(p) + distrib(gen));
        r2.emplace_back(p);
    }
    r2.emplace_back(0, 0);
    std::cout << "\n\n-----------------" << std::endl;
    auto before = std::chrono::system_clock::now();
    std::cout << "run self r1 ----------------------:" << std::endl;
    auto sr1 = self_or(r1);
    std::cout << bg::wkt(sr1) << std::endl;
    std::cout << "run self r2 ----------------------:" << std::endl;
    auto sr2 = self_or(r2);
    std::cout << bg::wkt(sr2) << std::endl;
    std::cout << "run r1 + r2 ----------------------:" << std::endl;
    auto sum = add(sr1, sr2);
    std::cout << bg::wkt(sum) << std::endl;
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}


void test_union(std::string first_s, std::string second_s, std::string ret_s) {
    multi_polygon first, second, ret;
    bg::read_wkt(first_s, first);
    assert(bg::is_valid(first));
    bg::read_wkt(second_s, second);
    assert(bg::is_valid(second));
    bg::read_wkt(ret_s, ret);
    assert(bg::is_valid(ret));
    std::cout << bg::wkt(add(first, second)) << std::endl;
    assert(bg::equals(add(first, second), ret));
}

void test_union_rectangle(int size) {
    using namespace std::chrono_literals;
    multi_polygon first, second;
    for (int i = 0; i < size; i++) {
        first.emplace_back(polygon{ {{0 + 2 * i, 0 + 2 * i}, {0 + 2 * i, 2 + 2 * i}, {2 + 2 * i, 2 + 2 * i}, {2 + 2 * i, 0 + 2 * i}, {0 + 2 * i, 0 + 2 * i}} });
        second.emplace_back(polygon{ {{1 + 2 * i, 1 + 2 * i}, {1 + 2 * i, 3 + 2 * i}, {3 + 2 * i, 3 + 2 * i}, {3 + 2 * i, 1 + 2 * i}, {1 + 2 * i, 1 + 2 * i}} });
    }
    auto before = std::chrono::system_clock::now();
    assert(bg::is_valid(first));
    assert(bg::is_valid(second));
    auto ret = add(first, second);
    assert(bg::is_valid(ret));
    assert(bg::area(ret) == 1 + 6 * size);
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

void test_self_or_rectangle(int size) {
    using namespace std::chrono_literals;
    multi_polygon poly;
    for (int i = 0; i < size; i++) {
        poly.emplace_back(polygon{ {{0 + i, 0 + i}, {0 + i, 2 + i}, {2 + i, 2 + i}, {2 + i, 0 + i}, {0 + i, 0 + i}} });
    }
    auto before = std::chrono::system_clock::now();
    std::cout << bg::wkt(poly) << std::endl;
    auto ret = self_or(poly);
    assert(bg::is_valid(ret));
    assert(bg::area(ret) == 1 + 3 * size);
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

int main()
{
    benchmark(400);
    /*
    benchmark(100);
    benchmark(200);
    benchmark(400);
    benchmark(1000);
    benchmark(2000);
    benchmark(4000);
    benchmark(10000);
    benchmark(20000);
    benchmark(40000);
    benchmark(80000);
    benchmark(100000);
    benchmark(200000);
    benchmark(400000);
    benchmark(1000000);
    benchmark(2000000);
    benchmark(4000000);
    benchmark(10000000);
    benchmark(20000000);
    benchmark(40000000);
    benchmark(100000000);
    */

    multi_polygon one, two, ret;

    test_union(
        "MULTIPOLYGON(((-59 867,-36 492,-182 486,-59 867)))",
        "MULTIPOLYGON(((-220 877,-54 821,-402 541,-808 638,-220 877)))",
        "MULTIPOLYGON(((-220 877,-72 827,-59 867,-56 822,-54 821,-56 819,-36 492,-182 486,-81 799,-402 541,-808 638,-220 877)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))",
        "MULTIPOLYGON(((2 2, 2 4, 4 4, 4 2, 2 2)))",
        "MULTIPOLYGON(((0 0, 0 3, 2 3, 2 4, 4 4, 4 2, 3 2, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 2, 5 1, 5 0, 0 0)))",
        "MULTIPOLYGON(((3 1, 0 2, 5 1, 3 1)))",
        "MULTIPOLYGON(((0 2,3 1,5 1,5 0,0 0,0 2)))"
    );
    test_union(
        "MULTIPOLYGON(((-1 -1, -1 3, 3 3, 3 -1, -1 -1), (0 0, 2 0, 2 2, 0 2, 0 0)))",
        "MULTIPOLYGON(((1 1, 1 4, 4 4, 4 1, 1 1)))",
        "MULTIPOLYGON(((-1 3,1 3,1 4,4 4,4 1,3 1,3 -1,-1 -1,-1 3),(2 0,2 1,1 1,1 2,0 2,0 0,2 0)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1), (6 6, 8 6, 8 8, 6 8, 6 6)))",
        "MULTIPOLYGON(((2 2, 2 7, 7 7, 7 2, 2 2)))",
        "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 2, 2 2, 2 3, 1 3, 1 1), (8 8, 6 8, 6 7, 7 7, 7 6, 8 6, 8 8)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 1 1, 2 1, 2 2, 3 3, 3 0, 0 0)))",
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 2 2, 1 2, 1 1, 0 0)))",
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))"
    );
    test_union(
        "MULTIPOLYGON(((-1461 -786,-1417 -833,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1030 -1246,-1608 -939,-1461 -786),(-1058 -1133,-1244 -963,-1301 -1039,-1058 -1133)),((-578 -670,-688 -679,-802 -443,-826 -400,-578 -670)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))",
        "MULTIPOLYGON(((-1450 -1280, -1450 -800, -1200 -1000, -1000 -1280, -1450 -1280)))",
        "MULTIPOLYGON(((-1461 -786,-1442 -807,-1410 -832,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1026 -1244,-1000 -1280,-1450 -1280,-1450 -1023,-1608 -939,-1461 -786),(-1228 -977,-1244 -963,-1245 -964,-1228 -977),(-1123 -1108,-1058 -1133,-1193 -1009,-1123 -1108)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))"
    );
        

    return 0;
}
