#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <ranges>
#include <optional>
#include <stack>
#include <random>
#include <chrono>
#include <boost/geometry.hpp>
#include <boost/container/flat_map.hpp>


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

auto count_time() {
    using namespace std::chrono_literals;
    static auto before_time = std::chrono::system_clock::now();
    auto after_time = std::chrono::system_clock::now();
    std::cout << "   before time: " << before_time << ", after time: " << after_time << ", runtime: " << (after_time - before_time) / 1s << "s" << std::endl;
    before_time = after_time;
}


auto construct_multi_polygons(auto&& rings) {
    std::vector<fake_bool> is_rings_cw(rings.size());
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (bg::area(rings[i]) > 0)
            is_rings_cw[i] = fake_bool::fake_true;
        else
            is_rings_cw[i] = fake_bool::fake_false;
    }
    std::vector<std::pair<box, std::size_t>> cw_rings_box;
    multi_polygon ret;

    cw_rings_box.reserve(rings.size());
    ret.reserve(rings.size());
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (is_rings_cw[i] == fake_bool::fake_true) {
            cw_rings_box.emplace_back(bg::return_envelope<box>(rings[i]), ret.size());
            polygon poly{ std::move(rings[i]) };
            ret.emplace_back(std::move(poly));
        }
    }
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
    // 判断 t1 和 t2 是否均在 [0, 1] 之间
    if (t1 > 0.0 && t1 < 1.0 && t2 > 0.0 && t2 < 1.0) {
        // -0.5 < x <= 0.5 to 0
        return point{ static_cast<ctype>(std::ceil(x1 + t1 * (x2 - x1) - 0.5)), static_cast<ctype>(std::ceil(y1 + t1 * (y2 - y1) - 0.5)) };
    }
    else {
        return {};
    }
}

// construct a graph with edge property (power) by segs
// segs should can create rings
auto construct_rings(const auto& segs, auto filter) {
    std::vector<std::pair<box, std::size_t>> boxes(segs.size());
    for (std::size_t i = 0; i < segs.size(); i++) {
        //std::cout << bg::wkt(segs[i]) << std::endl;
        boxes[i] = { bg::return_envelope<box>(segs[i]), i };
    }
    std::cout << "begin";
    count_time();
    bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > segs_box_rtree(boxes);

    std::vector<point> hot_pixels;
    hot_pixels.reserve(segs.size() * 2);

    for (const auto& seg : segs) {
        hot_pixels.emplace_back(bg::get<0, 0>(seg), bg::get<0, 1>(seg));
    }

    std::cout << "build segs rtree";
    count_time();
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
    std::cout << "find hot_pixels";
    count_time();
    {
        std::sort(std::begin(hot_pixels), std::end(hot_pixels),
            [](auto p1, auto p2) {
                // can dig more to see if we can use a better sort algorithm to better order pixels
                return std::pair{ bg::get<0>(p1), bg::get<1>(p1) } < std::pair{ bg::get<0>(p2), bg::get<1>(p2) };
            }
        ); // need define a better one
        auto last = std::unique(std::begin(hot_pixels), std::end(hot_pixels), [](auto p1, auto p2) { return bg::equals(p1, p2); });
        hot_pixels.erase(last, hot_pixels.end());
    }
    std::cout << "order hot pixels: ";
    count_time();

    std::vector<std::pair<std::size_t, std::size_t> > seg_pixel_pairs;
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        // have some precision problem, so expand 10
        constexpr auto expand = 1;
        auto min_corner = point{ bg::get<0>(hot_pixels[i]) - expand, bg::get<1>(hot_pixels[i]) - expand };
        auto max_corner = point{ bg::get<0>(hot_pixels[i]) + expand, bg::get<1>(hot_pixels[i]) + expand };

        std::for_each(segs_box_rtree.qbegin(bg::index::intersects(box{ min_corner, max_corner })), segs_box_rtree.qend(),
            [&](auto const& val) {
                // some more judeg, can be more precise
                auto x = bg::get<0>(hot_pixels[i]);
                auto y = bg::get<1>(hot_pixels[i]);
                auto x1 = bg::get<0, 0>(segs[val.second]);
                auto y1 = bg::get<0, 1>(segs[val.second]);
                auto x2 = bg::get<1, 0>(segs[val.second]);
                auto y2 = bg::get<1, 1>(segs[val.second]);
                if ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) < 0 || (x2 - x1) * (x2 - x) + (y2 - y1) * (y2 - y) < 0) return;

                if ((2 * x + 1 - x1 - x2) * (y1 - y2) > (2 * y + 1 - y1 - y2) * (x1 - x2) &&
                    (2 * x + 1 - x1 - x2) * (y1 - y2) >= (2 * y - 1 - y1 - y2) * (x1 - x2) &&
                    (2 * x - 1 - x1 - x2) * (y1 - y2) >= (2 * y + 1 - y1 - y2) * (x1 - x2) &&
                    (2 * x - 1 - x1 - x2) * (y1 - y2) >= (2 * y - 1 - y1 - y2) * (x1 - x2)) {
                    return;
                }
                if ((2 * x + 1 - x1 - x2) * (y1 - y2) < (2 * y + 1 - y1 - y2) * (x1 - x2) &&
                    (2 * x + 1 - x1 - x2) * (y1 - y2) <= (2 * y - 1 - y1 - y2) * (x1 - x2) &&
                    (2 * x - 1 - x1 - x2) * (y1 - y2) <= (2 * y + 1 - y1 - y2) * (x1 - x2) &&
                    (2 * x - 1 - x1 - x2) * (y1 - y2) <= (2 * y - 1 - y1 - y2) * (x1 - x2)) {
                    return;
                }
                seg_pixel_pairs.emplace_back(val.second, i);
            }
        );

    }
    std::cout << "find segs on hot_pixels: ";
    count_time();

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
        std::sort(std::begin(seg_pixel_pairs), std::end(seg_pixel_pairs));
        auto cur_begin = std::begin(seg_pixel_pairs);
        auto cur_last = cur_begin;
        while (++cur_last != std::end(seg_pixel_pairs)) {
            if (std::next(cur_last) != std::end(seg_pixel_pairs) && cur_begin->first == std::next(cur_last)->first) continue;
            std::sort(cur_begin, std::next(cur_last), less_by_segment{ segs[cur_begin->first], hot_pixels }); // should sort by the order of point in seg

            for (; cur_begin != cur_last; cur_begin++) {
                edges.emplace_back(cur_begin->second, std::next(cur_begin)->second);
            }
            cur_begin = std::next(cur_last);
        }
    }
    std::cout << "build edges: ";
    count_time();

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

    std::cout << "build edges power";
    count_time();

    std::vector<std::pair<bool, std::size_t> > direct_edges(edges.size() * 2);
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

    for (std::size_t i = 0; i < edges.size(); i++) {
        direct_edges[i] = { true, i };
        direct_edges[edges.size() + i] = { false, i };
    }

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
    std::sort(std::begin(direct_edges), std::end(direct_edges),
        [&](auto i1, auto i2) {
            if (source(i1) != source(i2)) [[likely]]
                return source(i1) < source(i2);
            return get_direction(i1) < get_direction(i2);
            // compare pair is slow
            //return std::pair{ source(i1), get_direction(i1) } < std::pair{ source(i2), get_direction(i2) };
        }
    );
    std::cout << "sort direct edges";
    count_time();

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
    std::cout << "connect direct edges";
    count_time();

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

    std::cout << "build faces";
    count_time();

    boost::container::flat_multimap<std::size_t, std::size_t> face_contain_relations;
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
            face_contain_relations.insert({ 0, face_id });
        else {
            for (; itr != faces_rtree.qend(); ++itr) {
                if (bg::relate(p, cw_faces[itr->second].first, bg::de9im::static_mask<'T', 'F', 'F'>{})) {
                    face_contain_relations.insert({ cw_faces[itr->second].second , face_id });
                    break;
                }
            }
        }
    }
    std::cout << "build face contain relations: ";
    count_time();

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
    std::cout << "traversal faces";
    count_time();

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

    // reomve co-point
    auto swap_co_source = [&](auto de1, auto de2) {
        auto pre1 = pre_direct_edges[de1];
        auto pre2 = pre_direct_edges[de2];
        next_direct_edges[pre1] = de2;
        next_direct_edges[pre2] = de1;
        pre_direct_edges[de1] = pre2;
        pre_direct_edges[de2] = pre1;
        };
    
    {
        std::vector<fake_bool> pixels_used(hot_pixels.size());
        for (std::size_t i = 0; i < direct_edges_exist.size(); i++) {
            if (direct_edges_exist[i] == fake_bool::fake_false) continue;
            auto cur_first_id = i;
            do {
                if (pixels_used[source(direct_edges[i])] == fake_bool::fake_false) {
                    pixels_used[source(direct_edges[i])] = fake_bool::fake_true;
                    i = pre_direct_edges[i];
                    continue;
                }
                else {
                    auto cur_begin_id = i;
                    auto begin_pixel_id = source(direct_edges[cur_begin_id]);
                    do {
                        i = next_direct_edges[i];
                        pixels_used[source(direct_edges[i])] = fake_bool::fake_false;
                    } while (source(direct_edges[i]) != source(direct_edges[cur_begin_id]));
                    swap_co_source(cur_begin_id, i);
                }
            } while (cur_first_id != i);
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
auto self_or(ring r) {
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

auto test(int size) {
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
    //std::cout << bg::wkt(sr1) << std::endl;
    std::cout << "run self r2 ----------------------:" << std::endl;
    self_or(r2);
    auto sr2 = self_or(r2);
    //std::cout << bg::wkt(sr2) << std::endl;
    std::cout << "run r1 + r2 ----------------------:" << std::endl;
    auto sum = add(r1, r2);
    //std::cout << bg::wkt(sum) << std::endl;
    auto after = std::chrono::system_clock::now();
    std::cout << "size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

int main()
{
    
    test(100);
    test(200);
    test(400);
    test(1000);
    test(2000);
    test(4000);
    test(10000);
    test(20000);
    test(40000);
    test(80000);
    test(100000);
    test(200000);
    test(400000);
    test(1000000);
    test(2000000);
    test(4000000);
    test(10000000);
    /*
    test(20000000);
    test(40000000);
    test(100000000);
    */

    multi_polygon one, two;

    bg::read_wkt(
        "MULTIPOLYGON(((0 0, 0 2, 5 1, 5 0, 0 0)))", one);
    assert(bg::is_valid(one));

    bg::read_wkt(
        "MULTIPOLYGON(((3 1, 0 2, 5 1, 3 1)))", two); // , ((0 0, 0 2, 5 0, 0 0))
    assert(bg::is_valid(two));
    //add(one, two);

    multi_polygon _one, _two;

    bg::read_wkt(
        "MULTIPOLYGON(((-1 -1, -1 3, 3 3, 3 -1, -1 -1), (0 0, 2 0, 2 2, 0 2, 0 0)))", _one);
    assert(bg::is_valid(_one));

    bg::read_wkt(
        "MULTIPOLYGON(((1 1, 1 4, 4 4, 4 1, 1 1)))", _two);
    assert(bg::is_valid(_two));

    add(_one, _two);

    return 0;
}
