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


    std::vector<std::tuple<std::size_t, std::size_t> > edges;
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
                edges.emplace_back(*cur_begin, *std::next(cur_begin));
            }
        }
    }
    log_done_time("build edges");
    return std::tuple{ std::move(edges), std::move(hot_pixels) };
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

struct duplicated_edge_t {
    auto source(const auto& edges) {
        if (i % 2 == 0) {
            return std::get<1>(edges[i / 2]);
        }
        else {
            return std::get<0>(edges[i / 2]);
        }
    }
    auto target(const auto& edges) {
        if (i % 2 == 0) {
            return std::get<0>(edges[i / 2]);
        }
        else {
            return std::get<1>(edges[i / 2]);
        }
    }
    template <auto power_index>
    auto power(const auto& edges) {
        if (i % 2 == 0) {
            return -std::get<power_index>(edges[i / 2]);
        }
        else {
            return std::get<power_index>(edges[i / 2]);
        }
    }
    auto dual() {
        if (i % 2 == 0) {
            return i + 1;
        }
        else {
            return i - 1;
        }
    }
    operator const std::size_t&() const{
        return i;
    }

    std::size_t i;
};

auto connect_duplicated_edges(const auto& edges, const auto& hot_pixels) {
    std::vector<duplicated_edge_t> duplicated_edges(edges.size() * 2);
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
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        auto cur_begin = std::begin(sort_duplicated_edges) + begin_location[i];
        auto cur_end = std::begin(sort_duplicated_edges) + end_location[i];
        if (cur_begin == cur_end) continue;
        assert(cur_begin + 1 != cur_end);
        next_duplicated_edges[(cur_end - 1)->dual()] = *cur_begin;
        for (++cur_begin; cur_begin != cur_end; ++cur_begin) {
            next_duplicated_edges[(cur_begin - 1)->dual()] = *cur_begin;
        }
    }
    return next_duplicated_edges;
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
            face_traversal.begin_edge({i});
            duplicated_edges_visited[i] = true;
            i = next_duplicated_edges[i];
        } while (i != cur_first);
    }
}

auto build_face_nearby_relations(const auto& next_duplicated_edges, const auto& edges_with_power) {
    struct face_id_calculation_traversal {
        void begin_face(std::size_t face_id) { cur_face_id = face_id; }
        void begin_edge(duplicated_edge_t duplicated_edge) {
            duplicated_edges_face_id[duplicated_edge] = cur_face_id;
        }
        std::vector<std::size_t>& duplicated_edges_face_id;
        std::size_t cur_face_id;
    };
    std::vector<std::size_t> duplicated_edges_face_id(next_duplicated_edges.size());
    traversal_face(next_duplicated_edges, face_id_calculation_traversal{ duplicated_edges_face_id, 0 });
    
    std::vector<std::tuple<std::size_t, std::size_t, int> > face_nearby_relations(next_duplicated_edges.size());
    for (std::size_t i = 0; i < face_nearby_relations.size(); i++) {
        duplicated_edge_t de{ i };
        face_nearby_relations[de] = { duplicated_edges_face_id[de.dual()], duplicated_edges_face_id[de], de.power<2>(edges_with_power)};
    }
}

// construct a graph with edge property (power) by segs
// segs should can create rings
auto construct_rings(auto segs, auto filter) {
    auto [edges, hot_pixels] = construct_graph(std::move(segs));
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
    {
#ifndef NDEBUG
        // check whether construct graph result is right
        // current geometry function is not precise, may have problems, here.
        std::vector<int> hot_pixels_times(hot_pixels.size());
        std::vector<int> hot_pixels_power(hot_pixels.size());
        for (auto edge_with_power : edges_with_power) {
            auto s = std::get<0>(edge_with_power);
            auto t = std::get<1>(edge_with_power);
            auto p = std::get<2>(edge_with_power);
            hot_pixels_times[s]++;
            hot_pixels_times[t]++;
            hot_pixels_power[s]+=p;
            hot_pixels_power[t]-=p;
        }
        for(std::size_t i = 0; i <hot_pixels.size(); i++) {
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
    // will use later
    //build_face_nearby_relations(connect_duplicated_edges(edges_with_power, hot_pixels), edges_with_power);

    auto source = [&](auto de) {
        if (de.first) return std::get<0>(edges_with_power[de.second]);
        else return std::get<1>(edges_with_power[de.second]);
        };
    auto target = [&](auto de) {
        if (de.first) return std::get<1>(edges_with_power[de.second]);
        else return std::get<0>(edges_with_power[de.second]);
        };
    auto power = [&](auto de) {
        if (de.first) return std::get<2>(edges_with_power[de.second]);
        else return -std::get<2>(edges_with_power[de.second]);
        };

    std::vector<std::pair<bool, std::size_t> > duplicated_edges(edges_with_power.size() * 2);
    {
        for (std::size_t i = 0; i < edges_with_power.size(); i++) {
            duplicated_edges[2 * i] = { true, i };
            duplicated_edges[2 * i + 1] = { false, i };
        }
        auto [begin_location, end_location, _duplicated_edges] = bucket_sort(
            std::move(duplicated_edges),
            hot_pixels.size(),
            [&](auto de) { return source(de); },
            [&](auto de) { return de; }
        );
        duplicated_edges = std::move(_duplicated_edges);
        for (std::size_t i = 0; i < hot_pixels.size(); i++) {
            if (end_location[i] - begin_location[i] < 3) continue;
            auto begin = std::begin(duplicated_edges);
            std::sort(begin + begin_location[i], begin + end_location[i],
                [&](auto i1, auto i2) {
                    return less_by_direction(hot_pixels[i], hot_pixels[target(i1)], hot_pixels[target(i2)] );
                }
            );
        }
    }
    log_done_time("sort direct edges");

    std::vector<std::pair<std::size_t, std::size_t> > direct_edge_pairs(edges_with_power.size());
    for (std::size_t i = 0; i < duplicated_edges.size(); i++) {
        if (duplicated_edges[i].first) direct_edge_pairs[duplicated_edges[i].second].second = i;
        else direct_edge_pairs[duplicated_edges[i].second].first = i;
    }
    auto get_dual = [&](auto i) {
        if (duplicated_edges[i].first) {
            return direct_edge_pairs[duplicated_edges[i].second].first;
        }
        else {
            return direct_edge_pairs[duplicated_edges[i].second].second;
        }
        };

    for (std::size_t i = 0; i < duplicated_edges.size(); i++) {
        auto dual = get_dual(i);
        //std::cout << "i = " << i << ", dual i = " << dual << std::endl;
    }

    // link direct edges, except the outer face, all faces are cw oriented
    std::vector<std::size_t> next_direct_edges(duplicated_edges.size());
    std::vector<std::size_t> pre_direct_edges(duplicated_edges.size());
    {
        auto cur_begin = std::begin(duplicated_edges);
        while (cur_begin != std::end(duplicated_edges)) {
            auto cur_end = not_adjacent_find(cur_begin, std::end(duplicated_edges), [&](auto de1, auto de2) { return source(de1) == source(de2); });
            pre_direct_edges[cur_begin - std::begin(duplicated_edges)] = get_dual(cur_end - 1 - std::begin(duplicated_edges));
            next_direct_edges[get_dual(cur_end - 1 - std::begin(duplicated_edges))] = cur_begin - std::begin(duplicated_edges);
            for (++cur_begin; cur_begin != cur_end; ++cur_begin) {
                pre_direct_edges[cur_begin - std::begin(duplicated_edges)] = get_dual(cur_begin - 1 - std::begin(duplicated_edges));
                next_direct_edges[get_dual(cur_begin - 1 - std::begin(duplicated_edges))] = cur_begin - std::begin(duplicated_edges);
            }
        }
    }
    log_done_time("connect direct edges");

    std::vector<std::size_t> edges_face_id(duplicated_edges.size());
    std::vector<std::pair<ring, std::size_t> > cw_faces;
    std::size_t cur_face_id = 1;
    boost::container::flat_map<std::size_t, std::size_t> ccw_face_id_to_direct_edge_id;
    for (std::size_t i = 0; i < duplicated_edges.size(); i++) {
        if (edges_face_id[i] != 0) continue;
        auto cur_first_id = i;
        ring r;
        // ring is closed, need to push one duplicated point
        r.push_back(hot_pixels[target(duplicated_edges[cur_first_id])]);
        do {
            edges_face_id[i] = cur_face_id;
            i = next_direct_edges[i];
            r.push_back(hot_pixels[target(duplicated_edges[i])]);
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
        std::vector<std::pair<std::size_t, std::size_t> > _face_contain_relations;
        std::vector<std::pair<box, std::size_t> > cw_faces_box(cw_faces.size());
        for (std::size_t i = 0; i < cw_faces.size(); i++) {
            cw_faces_box[i].first = bg::return_envelope<box>(cw_faces[i].first);
            cw_faces_box[i].second = i;
        }
        bg::index::rtree< std::pair<box, std::size_t>, bg::index::quadratic<128> > faces_rtree{ cw_faces_box };
        for (auto [face_id, direct_edge_id] : ccw_face_id_to_direct_edge_id) {
            point p = hot_pixels[source(duplicated_edges[direct_edge_id])];
            auto itr = faces_rtree.qbegin(bg::index::contains(p));
            if (itr == faces_rtree.qend())
                _face_contain_relations.emplace_back(0, face_id);
            else {
                for (; itr != faces_rtree.qend(); ++itr) {
                    if (bg::relate(p, cw_faces[itr->second].first, bg::de9im::static_mask<'T', 'F', 'F'>{})) {
                        _face_contain_relations.emplace_back(cw_faces[itr->second].second, face_id);
                        break;
                    }
                }
            }
        }
        face_contain_relations.insert(std::begin(_face_contain_relations), std::end(_face_contain_relations));
    }
    cw_faces.clear();
    log_done_time("build face contain relations");

    std::vector<bool> direct_edges_exist(duplicated_edges.size());
    {
        std::vector<std::optional<int> > faces_cw_power(cur_face_id);
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
                if (filter(faces_cw_power[cur_face_id].value())) direct_edges_exist[cur_direct_edge] = true;
                else direct_edges_exist[cur_direct_edge] = false;

                auto dual = get_dual(cur_direct_edge);
                auto next_face_id = edges_face_id[dual];
                if (!faces_cw_power[next_face_id]) {
                    faces_cw_power[next_face_id] = faces_cw_power[cur_face_id].value() + power(duplicated_edges[dual]);
                    stk.push({ next_face_id, dual });
                }
                cur_direct_edge = next_direct_edges[cur_direct_edge];
            } while (cur_direct_edge != cur_first_direct_edge);
        }
    }
    ccw_face_id_to_direct_edge_id.clear();
    log_done_time("traversal faces");

    for (auto&& [de1, de2] : direct_edge_pairs) {
        if (direct_edges_exist[de1] == true && direct_edges_exist[de2] == true) {
            direct_edges_exist[de1] = false;
            direct_edges_exist[de2] = false;
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

            if (direct_edges_exist[i] == false) continue;

            std::size_t cur_first_id = i;
            ring r;
            // ring is closed, need to push one duplicated point
            r.push_back(hot_pixels[target(duplicated_edges[i])]);
            direct_edges_exist[i] = false;
            do {
                i = next_direct_edges[i];
                r.push_back(hot_pixels[target(duplicated_edges[i])]);
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

auto benchmark(int size) {
    using namespace std::chrono_literals;
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(-600, 600);

    ring r1, r2;
    point p{ 0, 0 };
    r1.emplace_back(0, 0);
    for (int i = 0; i < size; i++) {
        bg::set<0>(p, bg::get<0>(p) + distrib(gen));
        bg::set<1>(p, bg::get<1>(p) + distrib(gen));
        r1.emplace_back(p);
    }
    r1.emplace_back(0, 0);
    p = point{ 0, 0 };
    r2.emplace_back(0, 0);
    for (int i = 0; i < size; i++) {
        bg::set<0>(p, bg::get<0>(p) + distrib(gen));
        bg::set<1>(p, bg::get<1>(p) + distrib(gen));
        r2.emplace_back(p);
    }
    r2.emplace_back(0, 0);
    std::cout << "\n\n-----------------" << std::endl;
    auto before = std::chrono::system_clock::now();
    std::cout << "run self r1 ----------------------:" << std::endl;
    std::cout << bg::wkt(r1) << std::endl;
    auto sr1 = self_or(r1);
    std::cout << bg::wkt(sr1) << std::endl;
    std::cout << "run self r2 ----------------------:" << std::endl;
    std::cout << bg::wkt(r2) << std::endl;
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
    assert(bg::area(ret) == (1 + 6 * size));
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
    assert(bg::area(ret) == (1 + 3 * size));
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

int main()
{
    //while(1)
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
