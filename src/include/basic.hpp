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
#include <type_traits>
#include <boost/geometry.hpp>
#include <boost/container/flat_map.hpp>
extern "C" {
#define _CRT_USE_C_COMPLEX_H
#include <GraphBLAS.h>
#undef _CRT_USE_C_COMPLEX_H
}

//#define LOG_TIME

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<int>;
using segment = bg::model::segment<point>;
using box = bg::model::box<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;

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
    std::vector<bool> is_rings_cw(rings.size());
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (bg::area(rings[i]) > 0)
            is_rings_cw[i] = true;
        else
            is_rings_cw[i] = false;
    }
    std::vector<std::pair<box, std::size_t>> cw_rings_box;
    multi_polygon ret;
    log_done_time("judge cw");

    cw_rings_box.reserve(rings.size());
    ret.reserve(rings.size());
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (is_rings_cw[i] == true) {
            cw_rings_box.emplace_back(bg::return_envelope<box>(rings[i]), ret.size());
            polygon poly{ std::move(rings[i]) };
            ret.emplace_back(std::move(poly));
        }
    }
    log_done_time("build multi polygon");

    bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > rings_box_tree(cw_rings_box);
    for (std::size_t i = 0; i < rings.size(); i++) {
        if (is_rings_cw[i] == false) {
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

struct less_by_segment {
    bool operator() (point p1, point p2) {
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
};

bool less_by_direction(point source, point target1, point target2) {
    constexpr auto get_direction = [](point v1, point v2) {
        // can be more precise
        double dx = bg::get<0>(v2) - bg::get<0>(v1);
        double dy = bg::get<1>(v2) - bg::get<1>(v1);
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
    return get_direction(source, target1) < get_direction(source, target2);
}

auto bucket_sort(auto vec, auto bucket_size, auto get_bucket, auto get_left) {
    std::vector<int> times(bucket_size);
    for (auto val : vec) {
        times[get_bucket(val)]++;
    }
    std::vector<unsigned int> begin_location(times.size());
    std::exclusive_scan(std::begin(times), std::end(times), std::begin(begin_location), 0);
    std::vector<std::invoke_result_t<decltype(get_left), typename decltype(vec)::value_type> > left(vec.size());
    auto current_location{ begin_location };
    for (auto val : vec) {
        left[current_location[get_bucket(val)]++] = get_left(val);
    }
    return std::tuple{ std::move(begin_location), std::move(current_location), std::move(left) };
}

// std has adjacent find, but the function result is not continuous
// it's better to use c++23 chunk_by
auto not_adjacent_find(auto begin, auto end, auto binary) {
    assert(begin != end);
    auto cur = begin;
    while (std::next(cur) != end) {
        if (binary(*cur, *std::next(cur))) {
            ++cur;
        }
        else {
            break;
        }
    }
    return std::next(cur);
}

auto construct_graph(auto segs) {
    log_done_time("start of construct_graph");
    std::vector<std::pair<box, std::size_t>> boxes(segs.size());
    for (std::size_t i = 0; i < segs.size(); i++) {
        boxes[i] = { bg::return_envelope<box>(segs[i]), i };
    }
    bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > segs_box_rtree(std::move(boxes));

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


    std::vector<std::size_t> edges_start;
    std::vector<std::size_t> edges_end;
    {
        auto [segs_begin_location, segs_end_location, pixels] = bucket_sort(
            seg_pixel_pairs,
            segs.size(),
            [](auto val) {return val.first; },
            [](auto val) {return val.second; }
        );
        for (std::size_t i = 0; i < segs.size(); i++) {
            auto cur_begin = std::begin(pixels) + segs_begin_location[i];
            auto cur_end = std::begin(pixels) + segs_end_location[i];
            auto cur_last = cur_end - 1;
            std::sort(cur_begin, cur_end,
                [&](auto pi, auto pj) {
                    return less_by_segment{ segs[i] }(hot_pixels[pi], hot_pixels[pj]);
                }
            );

            for (; cur_begin != cur_last; cur_begin++) {
                edges_start.emplace_back(*cur_begin);
                edges_end.emplace_back(* std::next(cur_begin));
            }
        }
    }
    log_done_time("build edges");
    return std::tuple{ std::move(edges_start), std::move(edges_end), std::move(hot_pixels) };
}

auto edges_direction_to_power(auto edges) {
    using edge_t = typename decltype(edges)::value_type;
    using new_edge_t = decltype(std::tuple_cat(std::declval<edge_t>(), std::tuple<int>{}));
    std::vector<new_edge_t> ret(edges.size());
    for (std::size_t i = 0; i < ret.size(); i++) {
        ret[i] = std::tuple_cat(edges[i], std::tuple{ 1 });
        auto& v1 = std::get<0>(ret[i]);
        auto& v2 = std::get<1>(ret[i]);
        if (v1 < v2) {
            std::get<std::tuple_size<new_edge_t>() - 1>(ret[i]) = 1;
        }
        else {
            std::swap(v1, v2);
            std::get<std::tuple_size<new_edge_t>() - 1>(ret[i]) = -1;
        }
    }
    return ret;
}

auto sort_edges(auto edges, auto vertex_number) {
    auto [begin_location, end_location, ordered_double_edges] = bucket_sort(
        std::move(edges),
        vertex_number,
        [](auto val) { return std::get<0>(val); },
        [](auto val) { return val; }
    );
    for (std::size_t i = 0; i < vertex_number; i++) {
        auto cur_begin = std::begin(ordered_double_edges) + begin_location[i];
        auto cur_end = std::begin(ordered_double_edges) + end_location[i];
        std::sort(cur_begin, cur_end, [](auto e1, auto e2) { return std::get<1>(e1) < std::get<1>(e2); });
    }
    return std::move(ordered_double_edges);
}

auto unique_edges(auto edges, auto merge_func) {
    constexpr auto equal = [](auto e1, auto e2) {
        return std::get<0>(e1) == std::get<0>(e2) && std::get<1>(e1) == std::get<1>(e2);
        };

    auto cur_begin = std::begin(edges);
    auto cur_result = cur_begin;
    while (cur_begin != std::end(edges)) {
        auto cur_end = not_adjacent_find(cur_begin, std::end(edges), equal);
        *cur_result++ = std::reduce(std::next(cur_begin), cur_end, *cur_begin, merge_func);
        cur_begin = cur_end;
    }
    edges.erase(cur_result, std::end(edges));
    return std::move(edges);
}

auto construct_edges_with_power(auto segs) {
    auto [_edges_start, _edges_end, hot_pixels] = construct_graph(std::move(segs));
	std::vector<std::pair<std::size_t, std::size_t> > edges;
	for (std::size_t i = 0; i < _edges_start.size(); i++) {
        edges.emplace_back(_edges_start[i], _edges_end[i]);
    }
    auto edges_with_power =
        unique_edges(
            sort_edges(edges_direction_to_power(std::move(edges)), hot_pixels.size()),
            [](auto e1, auto e2) {
                return std::tuple{ std::get<0>(e1), std::get<1>(e1), std::get<2>(e1) + std::get<2>(e2) };
            }
        );
    std::erase_if(edges_with_power, [](auto edge_with_power) {
        return (std::get<2>(edge_with_power) == 0);
        });
    std::vector<std::size_t> edges_start(edges_with_power.size());
    std::vector<std::size_t> edges_end(edges_with_power.size());
    std::vector<std::size_t> powers(edges_with_power.size());
	for (std::size_t i = 0; i < edges_with_power.size(); i++) {
        edges_start[i] = std::get<0>(edges_with_power[i]);
        edges_end[i] = std::get<1>(edges_with_power[i]);
        powers[i] = std::get<2>(edges_with_power[i]);
    }
    {
#ifndef NDEBUG
        // check whether construct graph result is right
        // current geometry function is not precise, may have problems, here.
        std::vector<int> hot_pixels_times(hot_pixels.size());
        std::vector<int> hot_pixels_power(hot_pixels.size());
        for (std::size_t i = 0; i < edges_start.size(); i++) {
            auto s = edges_start[i];
            auto t = edges_end[i];
            hot_pixels_times[s]++;
            hot_pixels_times[t]++;
            hot_pixels_power[s] += powers[i];
			hot_pixels_power[t] -= powers[i];
        }
        for (std::size_t i = 0; i < hot_pixels.size(); i++) {
            if (hot_pixels_times[i] == 1 || hot_pixels_power[i] != 0) {
                std::cout << hot_pixels_times[i] << std::endl;
                std::cout << hot_pixels_power[i] << std::endl;
                std::cout << bg::wkt(hot_pixels[i]) << std::endl;
                assert(false);
            }
        }
#endif
    }
    log_done_time("calculate edge power");
    return std::tuple{ std::move(hot_pixels), std::tuple{std::move(edges_start), std::move(edges_end), std::move(powers)} };
}

struct duplicated_edge_t {
    auto source(const auto& edges) {
        if (i % 2 == 0) {
            return std::get<1>(edges)[i / 2];
        }
        else {
            return std::get<0>(edges)[i / 2];
        }
    }
    auto target(const auto& edges) {
        if (i % 2 == 0) {
            return std::get<0>(edges)[i / 2];
        }
        else {
            return std::get<1>(edges)[i / 2];
        }
    }
    template <auto power_index>
    auto power(const auto& edges) {
        if (i % 2 == 0) {
            return -std::get<power_index>(edges)[i / 2];
        }
        else {
            return std::get<power_index>(edges)[i / 2];
        }
    }
    auto dual() {
        if (i % 2 == 0) {
            return duplicated_edge_t{ i + 1 };
        }
        else {
            return duplicated_edge_t{ i - 1 };
        }
    }
    operator const std::size_t& () const {
        return i;
    }

    std::size_t i;
};

auto connect_duplicated_edges(const auto& edges, const auto& hot_pixels) {
    std::vector<duplicated_edge_t> duplicated_edges(std::get<0>(edges).size() * 2);
    for (std::size_t i = 0; i < duplicated_edges.size(); i++) {
        duplicated_edges[i] = duplicated_edge_t{ i };
    }

    auto [begin_location, end_location, sort_duplicated_edges] = bucket_sort(
        std::move(duplicated_edges),
        hot_pixels.size(),
        [&](auto de) { return de.source(edges); },
        [&](auto de) { return de; }
    );
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        if (end_location[i] - begin_location[i] < 3) continue;
        std::sort(
            std::begin(sort_duplicated_edges) + begin_location[i],
            std::begin(sort_duplicated_edges) + end_location[i],
            [&](auto i1, auto i2) {
                return less_by_direction(hot_pixels[i], hot_pixels[i1.target(edges)], hot_pixels[i2.target(edges)]);
            }
        );
    }
    std::vector<duplicated_edge_t> next_duplicated_edges(sort_duplicated_edges.size());
    std::vector<duplicated_edge_t> pre_duplicated_edges(sort_duplicated_edges.size());
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        auto cur_begin = std::begin(sort_duplicated_edges) + begin_location[i];
        auto cur_end = std::begin(sort_duplicated_edges) + end_location[i];
        if (cur_begin == cur_end) continue;
        assert(cur_begin + 1 != cur_end);
        next_duplicated_edges[(cur_end - 1)->dual()] = *cur_begin;
        pre_duplicated_edges[*cur_begin] = (cur_end - 1)->dual();
        for (++cur_begin; cur_begin != cur_end; ++cur_begin) {
            next_duplicated_edges[(cur_begin - 1)->dual()] = *cur_begin;
            pre_duplicated_edges[*cur_begin] = (cur_begin - 1)->dual();
        }
    }
    return std::pair{ std::move(next_duplicated_edges), std::move(pre_duplicated_edges) };
}

auto traversal_face(const std::vector<duplicated_edge_t>& next_duplicated_edges, auto face_traversal) {
    auto size = next_duplicated_edges.size();
    std::vector<bool> duplicated_edges_visited(size);
    std::size_t cur_face_id = 0;

    for (std::size_t i = 0; i < size; i++) {
        if (duplicated_edges_visited[i]) continue;
        ++cur_face_id;
        duplicated_edge_t cur_first{ i };
        face_traversal.begin_face(cur_face_id);
        do {
            face_traversal.begin_edge({ i });
            duplicated_edges_visited[i] = true;
            i = next_duplicated_edges[i];
        } while (i != cur_first);
        face_traversal.end_face();
    }
}

auto log_duplicated_edges_face_id(const auto& next_duplicated_edges) {
    struct face_id_calculation_traversal {
        void begin_face(std::size_t face_id) { cur_face_id = face_id; }
        void begin_edge(duplicated_edge_t duplicated_edge) {
            duplicated_edges_face_id[duplicated_edge] = cur_face_id;
        }
        void end_face() {}
        std::vector<std::size_t>& duplicated_edges_face_id;
        std::size_t& cur_face_id;
    };

    std::vector<std::size_t> duplicated_edges_face_id(next_duplicated_edges.size());
    std::size_t cur_face_id = 1;
    face_id_calculation_traversal face_id_traversal{ duplicated_edges_face_id, cur_face_id };
    traversal_face(next_duplicated_edges, face_id_traversal);

    return std::pair{ std::move(duplicated_edges_face_id), cur_face_id + 1 };
}

auto build_face_nearby_relations(const auto& next_duplicated_edges, const auto& duplicated_edges_face_id) {
    std::vector<std::tuple<std::size_t, std::size_t, duplicated_edge_t> > face_nearby_relations(next_duplicated_edges.size());
    for (std::size_t i = 0; i < face_nearby_relations.size(); i++) {
        duplicated_edge_t de{ i };
        face_nearby_relations[de] = { duplicated_edges_face_id[de.dual()], duplicated_edges_face_id[de], de };
        face_nearby_relations[de] = { duplicated_edges_face_id[de], duplicated_edges_face_id[de.dual()], de.dual()};
    }

    log_done_time("build face nearby relations");

    return face_nearby_relations;
}

auto build_face_contains_relations(const auto& next_duplicated_edges, const auto& edges_with_power, const auto& hot_pixels) {
    struct faces_record_traversal {
        void begin_face(std::size_t face_id) {
            cur_face_id = face_id;
            cur_r = {};
        }
        void begin_edge(duplicated_edge_t duplicated_edge) {
            cur_r.push_back(_hot_pixels[duplicated_edge.target(_edges_with_power)]);
        }
        void end_face() {
            cur_r.push_back(cur_r[0]); // ring is closed, need to push one duplicated point
            if (bg::area(cur_r) > 0) { // for cw orientation, the area > 0, may use more precise algorithm
                cw_faces.emplace_back(std::move(cur_r), cur_face_id);
            }
            else {
                ccw_face_id_to_direct_edge.emplace_back(cur_r[0], cur_face_id);
            }
        }
        decltype(hot_pixels) _hot_pixels;
        decltype(edges_with_power) _edges_with_power;
        std::vector<std::pair<ring, std::size_t> >& cw_faces;
        std::vector<std::pair<point, std::size_t> >& ccw_face_id_to_direct_edge;
        std::size_t& cur_face_id;

        ring cur_r{};
    };

    std::vector<std::pair<ring, std::size_t> > cw_faces;
    std::vector<std::pair<point, std::size_t> > ccw_faces;
    std::size_t cur_face_id;
    traversal_face(next_duplicated_edges, faces_record_traversal{ hot_pixels, edges_with_power, cw_faces, ccw_faces, cur_face_id });
    auto face_num = cur_face_id + 1;

    std::vector<std::pair<std::size_t, std::size_t> > face_full_contain_relations;
    {
        std::vector<std::pair<box, std::size_t> > cw_faces_box(cw_faces.size());
        for (std::size_t i = 0; i < cw_faces.size(); i++) {
            cw_faces_box[i].first = bg::return_envelope<box>(cw_faces[i].first);
            cw_faces_box[i].second = i;
        }
        bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > faces_rtree{ cw_faces_box };
        for (auto [p, face_id] : ccw_faces) {
            auto itr = faces_rtree.qbegin(bg::index::contains(p));
            face_full_contain_relations.emplace_back(0, face_id);
            for (auto itr = faces_rtree.qbegin(bg::index::contains(p)); itr != faces_rtree.qend(); ++itr) {
                if (bg::relate(p, cw_faces[itr->second].first, bg::de9im::static_mask<'T', 'F', 'F'>{})) {
                    face_full_contain_relations.emplace_back(cw_faces[itr->second].second, face_id);
                }
            }
        }
    }
    std::vector<std::size_t> out_faces_times(face_num);
    for (auto [out_face, in_face] : face_full_contain_relations) {
        out_faces_times[out_face]++;
    }

    auto [begin_location, end_location, out_faces] = bucket_sort(
        std::move(face_full_contain_relations),
        face_num,
        [](auto pair) { return std::get<1>(pair); },
        [](auto pair) { return std::get<0>(pair); }
    );

    std::vector<std::pair<std::size_t, std::size_t> > face_contain_relations;

    for (std::size_t in_face = 0; in_face < face_num; in_face++) {
        if (end_location[in_face] - begin_location[in_face]  == 0) continue;
        auto out_face = out_faces[begin_location[in_face] ];
        for (std::size_t j = begin_location[in_face]; j < end_location[in_face]; j++) {
            if (out_faces_times[out_faces[j] ] < out_faces_times[out_face]) {
                out_face = out_faces[j];
            }
        }
        face_contain_relations.emplace_back(out_face, in_face);
    }

    log_done_time("build face contains relations");
    
    return face_contain_relations;
}

template <auto power_i>
auto group_face_relations(auto face_contain_relations, auto face_nearby_relations, auto face_num, const auto& edges_with_power) {
    std::vector<std::tuple<std::size_t, std::size_t, int> > face_relations(face_contain_relations.size() + face_nearby_relations.size());
    for (std::size_t i = 0; i < face_contain_relations.size(); i++) {
        auto t = std::tuple<std::size_t, std::size_t, int>{std::get<0>(face_contain_relations[i]), std::get<1>(face_contain_relations[i]), 0};
        face_relations[i] = t;
    }
    for (std::size_t i = 0; i < face_nearby_relations.size(); i++) {
        face_relations[i + face_contain_relations.size()] = std::tuple{
            std::get<0>(face_nearby_relations[i]),
            std::get<1>(face_nearby_relations[i]),
            std::get<2>(face_nearby_relations[i]).template power<power_i>(edges_with_power)
        };
    }
    return bucket_sort(
        face_relations,
        face_num,
        [](auto val) { return std::get<0>(val); },
        [](auto val) { return std::pair{ std::get<1>(val), std::get<2>(val) }; }
    );
    log_done_time("group face relations");
    
}

// construct a graph with edge property (power) by segs
// segs should can create rings
auto construct_rings(auto segs, auto filter) {
    auto [hot_pixels, edges_with_power] = construct_edges_with_power(std::move(segs));
    log_done_time("calculate edge power");

    auto [next_duplicated_edges, pre_duplicated_edges] = connect_duplicated_edges(edges_with_power, hot_pixels);
    log_done_time("connect direct edges");

    auto [duplicated_edges_face_id, face_num] = log_duplicated_edges_face_id(next_duplicated_edges);
    log_done_time("build faces");

    std::vector<bool> faces_exist(face_num);
    {
        auto [begin_location, end_location, face_relations] = group_face_relations<2>(
            build_face_contains_relations(next_duplicated_edges, edges_with_power, hot_pixels),
            build_face_nearby_relations(next_duplicated_edges, duplicated_edges_face_id),
            face_num,
            edges_with_power
        );

        std::vector<int > faces_cw_power(face_num);
        std::vector<bool> faces_visited(face_num);
        std::stack<std::size_t > stk;
        faces_visited[0] = true;
        faces_cw_power[0] = 0;
        stk.push(0);
        while (!stk.empty()) {
            auto cur_face_id = stk.top();
            stk.pop();

            if (filter(faces_cw_power[cur_face_id]) ) faces_exist[cur_face_id] = true;

            for (auto i = begin_location[cur_face_id]; i < end_location[cur_face_id]; i++) {
                auto [next_face_id, power] = face_relations[i];
                if (!faces_visited[next_face_id]) {
                    faces_visited[next_face_id] = true;
                    faces_cw_power[next_face_id] = faces_cw_power[cur_face_id] + power;
                    stk.push({ next_face_id });
                }
            }
        }
    }
    log_done_time("traversal faces");

    std::vector<bool> direct_edges_exist(std::get<0>(edges_with_power).size() * 2);
    for (std::size_t i = 0; i < std::get<0>(edges_with_power).size() * 2; i++) {
        direct_edges_exist[i] = faces_exist[duplicated_edges_face_id[i]];
    }
    for (std::size_t i = 0; i < std::get<0>(edges_with_power).size(); i++) {
        auto de1 = duplicated_edge_t{ 2 * i };
        auto de2 = duplicated_edge_t{ 2 * i + 1 };
        if (direct_edges_exist[de1] == true && direct_edges_exist[de2] == true) {
            direct_edges_exist[de1] = false;
            direct_edges_exist[de2] = false;
            auto pre1 = pre_duplicated_edges[de1];
            auto pre2 = pre_duplicated_edges[de2];
            auto next1 = next_duplicated_edges[de1];
            auto next2 = next_duplicated_edges[de2];
            pre_duplicated_edges[next1] = pre2;
            pre_duplicated_edges[next2] = pre1;
            next_duplicated_edges[pre1] = next2;
            next_duplicated_edges[pre2] = next1;
            // no use in fact
            pre_duplicated_edges[de1] = de2;
            pre_duplicated_edges[de2] = de1;
            next_duplicated_edges[de1] = de2;
            next_duplicated_edges[de2] = de1;
        }
    }

    std::vector<ring> ret_rings;
    {

        for (std::size_t i = 0; i < direct_edges_exist.size(); i++) {

            if (direct_edges_exist[i] == false) continue;

            std::size_t cur_first_id = i;
            ring r;
            // ring is closed, need to push one duplicated point
            r.push_back(hot_pixels[duplicated_edge_t{ i }.target(edges_with_power)]);
            direct_edges_exist[i] = false;
            do {
                i = next_duplicated_edges[i];
                r.push_back(hot_pixels[duplicated_edge_t{ i }.target(edges_with_power)]);
                direct_edges_exist[i] = false;
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

auto add(const auto& ps1, const auto& ps2) {
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
    return construct_multi_polygons(construct_rings(std::move(segs), filter));
}

auto intersection(const auto& ps1, const auto& ps2) {
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

    constexpr auto filter = [](auto cw_power) { return cw_power > 1; };
    return construct_multi_polygons(construct_rings(std::move(segs), filter));
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
    auto rings = construct_rings(std::move(segs), [](auto cw_power) { return cw_power > 0; });
    auto ret = construct_multi_polygons(std::move(rings));
    return ret;
}