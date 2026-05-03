#pragma once

#include <algorithm>
#include <boost/container/flat_map.hpp>
#include <boost/geometry.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <vector>

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<int>;
using segment = bg::model::segment<point>;
using box = bg::model::box<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

inline auto bucket_sort(auto vec, auto bucket_size, auto get_bucket, auto get_left) {
    std::vector<int> times(bucket_size);
    for (auto val : vec) times[get_bucket(val)]++;
    std::vector<unsigned int> begin_location(times.size());
    std::exclusive_scan(std::begin(times), std::end(times), std::begin(begin_location), 0);
    std::vector<std::invoke_result_t<decltype(get_left), typename decltype(vec)::value_type>> left(
        vec.size());
    auto current_location{begin_location};
    for (auto val : vec) left[current_location[get_bucket(val)]++] = get_left(val);
    return std::tuple{std::move(begin_location), std::move(current_location), std::move(left)};
}

inline auto not_adjacent_find(auto begin, auto end, auto binary) {
    assert(begin != end);
    auto cur = begin;
    while (std::next(cur) != end) {
        if (binary(*cur, *std::next(cur))) ++cur;
        else break;
    }
    return std::next(cur);
}

// Undirected graph connected components via DFS. O(n + m).
inline std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>>& edges) {

    std::vector<std::pair<std::size_t, std::size_t>> all_edges;
    all_edges.reserve(edges.size() * 2);
    for (auto [a, b] : edges) {
        all_edges.emplace_back(a, b);
        all_edges.emplace_back(b, a);
    }

    auto [adj_begin, adj_end, adjacency] = bucket_sort(
        all_edges, n,
        [](const std::pair<std::size_t, std::size_t>& e) { return e.first; },
        [](const std::pair<std::size_t, std::size_t>& e) { return e.second; });

    std::vector<std::size_t> comp(n, ~0ULL);
    std::vector<std::size_t> stack;
    std::size_t num_components = 0;

    for (std::size_t v = 0; v < n; v++) {
        if (comp[v] != ~0ULL) continue;
        comp[v] = num_components++;
        stack.push_back(v);
        while (!stack.empty()) {
            auto u = stack.back(); stack.pop_back();
            for (auto j = adj_begin[u]; j < adj_end[u]; j++) {
                auto w = adjacency[j];
                if (comp[w] == ~0ULL) {
                    comp[w] = comp[v];
                    stack.push_back(w);
                }
            }
        }
    }

    return comp;
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

struct less_by_segment {
    bool operator()(point p1, point p2) const {
        int64_t x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
        int64_t x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
        int64_t x3 = bg::get<0, 0>(s), y3 = bg::get<0, 1>(s);
        int64_t x4 = bg::get<1, 0>(s), y4 = bg::get<1, 1>(s);
        return (x2 - x1) * (x4 - x3) + (y2 - y1) * (y4 - y3) > 0;
    }
    segment s;
};

inline bool less_by_direction(point source, point target1, point target2) {
    auto get_direction = [](point v1, point v2) {
        int64_t dx = bg::get<0>(v2) - bg::get<0>(v1);
        int64_t dy = bg::get<1>(v2) - bg::get<1>(v1);
        enum class quadrant { _1, _2, _3, _4, zero };
        if (dx > 0 && dy >= 0) return std::pair{quadrant::_1, (long double)dy / dx};
        if (dx <= 0 && dy > 0) return std::pair{quadrant::_2, -(long double)dx / dy};
        if (dx < 0 && dy <= 0) return std::pair{quadrant::_3, (long double)dy / dx};
        if (dx >= 0 && dy < 0) return std::pair{quadrant::_4, -(long double)dx / dy};
        return std::pair{quadrant::zero, 0.0L};
    };
    return get_direction(source, target1) < get_direction(source, target2);
}

inline int64_t cross(int64_t x1, int64_t y1, int64_t x2, int64_t y2) {
    return x1 * y2 - y1 * x2;
}

inline std::optional<point> get_intersection(segment s1, segment s2) {
    int64_t x1 = bg::get<0, 0>(s1), y1 = bg::get<0, 1>(s1);
    int64_t x2 = bg::get<1, 0>(s1), y2 = bg::get<1, 1>(s1);
    int64_t x3 = bg::get<0, 0>(s2), y3 = bg::get<0, 1>(s2);
    int64_t x4 = bg::get<1, 0>(s2), y4 = bg::get<1, 1>(s2);
    int64_t dx1 = x2 - x1, dy1 = y2 - y1;
    int64_t dx2 = x4 - x3, dy2 = y4 - y3;
    int64_t d = cross(dx1, dy1, dx2, dy2);
    if (d == 0) return {};
    int64_t t1_num = cross(x3 - x1, y3 - y1, dx2, dy2);
    int64_t t2_num = cross(x3 - x1, y3 - y1, dx1, dy1);
    if (d < 0) { d = -d; t1_num = -t1_num; t2_num = -t2_num; }
    if (t1_num <= 0 || t1_num >= d || t2_num <= 0 || t2_num >= d) return {};
    long double fx = (long double)x1 + (long double)t1_num * dx1 / d;
    long double fy = (long double)y1 + (long double)t1_num * dy1 / d;
    return point{static_cast<int>(std::ceil(fx - 0.5L)),
                 static_cast<int>(std::ceil(fy - 0.5L))};
}

inline bool is_point_on_segment(point p, segment s) {
    int64_t x = bg::get<0>(p), y = bg::get<1>(p);
    int64_t x1 = bg::get<0, 0>(s), y1 = bg::get<0, 1>(s);
    int64_t x2 = bg::get<1, 0>(s), y2 = bg::get<1, 1>(s);
    if ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) < 0) return false;
    if ((x2 - x1) * (x2 - x) + (y2 - y1) * (y2 - y) < 0) return false;
    int64_t dx = x2 - x1, dy = y2 - y1;
    if ((2 * x + 1 - x1 - x2) * dy > (2 * y + 1 - y1 - y2) * dx &&
        (2 * x + 1 - x1 - x2) * dy >= (2 * y - 1 - y1 - y2) * dx &&
        (2 * x - 1 - x1 - x2) * dy >= (2 * y + 1 - y1 - y2) * dx &&
        (2 * x - 1 - x1 - x2) * dy >= (2 * y - 1 - y1 - y2) * dx)
        return false;
    if ((2 * x + 1 - x1 - x2) * dy < (2 * y + 1 - y1 - y2) * dx &&
        (2 * x + 1 - x1 - x2) * dy <= (2 * y - 1 - y1 - y2) * dx &&
        (2 * x - 1 - x1 - x2) * dy <= (2 * y + 1 - y1 - y2) * dx &&
        (2 * x - 1 - x1 - x2) * dy <= (2 * y - 1 - y1 - y2) * dx)
        return false;
    return true;
}

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

struct edge_t { std::size_t start; std::size_t end; };
struct edge_with_power_t { std::size_t start, end; int power; };

struct chain_build_result {
    std::vector<std::size_t> indices, offsets;
    std::vector<int> powers;
    std::vector<std::size_t> edge_to_chain;
};

struct half_chain {
    std::size_t id;
    std::size_t chain_id() const { return id / 2; }
    bool is_forward() const { return id % 2 == 0; }
    std::size_t source_node(const chain_build_result& c) const {
        if (is_forward()) return c.indices[c.offsets[chain_id()]];
        else return c.indices[c.offsets[chain_id() + 1] - 1];
    }
    std::size_t target_node(const chain_build_result& c) const {
        if (is_forward()) return c.indices[c.offsets[chain_id() + 1] - 1];
        else return c.indices[c.offsets[chain_id()]];
    }
    std::size_t next_along_source(const chain_build_result& c) const {
        if (is_forward()) return c.indices[c.offsets[chain_id()] + 1];
        else return c.indices[c.offsets[chain_id() + 1] - 2];
    }
    int power(const chain_build_result& c) const {
        if (is_forward()) return c.powers[chain_id()];
        else return -c.powers[chain_id()];
    }
    half_chain dual() const { return {id ^ 1}; }
};

// ---------------------------------------------------------------------------
// Step 1-3: edges with power
// ---------------------------------------------------------------------------

inline std::tuple<std::vector<edge_t>, std::vector<point>> construct_graph(auto segments) {
    std::vector<std::pair<box, std::size_t>> boxes(segments.size());
    for (std::size_t i = 0; i < segments.size(); i++)
        boxes[i] = {bg::return_envelope<box>(segments[i]), i};
    bg::index::rtree<std::pair<box, std::size_t>, bg::index::quadratic<128>> segments_box_rtree(
        std::move(boxes));

    std::vector<point> hot_pixels;
    hot_pixels.reserve(segments.size() * 2);
    for (const auto& segment : segments) {
        hot_pixels.emplace_back(bg::get<0, 0>(segment), bg::get<0, 1>(segment));
        hot_pixels.emplace_back(bg::get<1, 0>(segment), bg::get<1, 1>(segment));
    }

    for (std::size_t i = 0; i < segments.size(); i++)
        std::for_each(
            segments_box_rtree.qbegin(bg::index::intersects(boxes[i].first)), segments_box_rtree.qend(),
            [&](auto const& o) {
                if (auto p = get_intersection(segments[i], segments[o.second]))
                    hot_pixels.push_back(p.value());
            });

    std::sort(std::begin(hot_pixels), std::end(hot_pixels), [](auto p1, auto p2) {
        return std::pair{bg::get<0>(p1), bg::get<1>(p1)} <
               std::pair{bg::get<0>(p2), bg::get<1>(p2)};
    });
    auto last = std::unique(std::begin(hot_pixels), std::end(hot_pixels),
                            [](auto p1, auto p2) { return bg::equals(p1, p2); });
    hot_pixels.erase(last, hot_pixels.end());

    std::vector<std::pair<std::size_t, std::size_t>> segment_pixel_pairs;
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        constexpr auto expand = 1;
        auto min_corner = point{bg::get<0>(hot_pixels[i]) - expand, bg::get<1>(hot_pixels[i]) - expand};
        auto max_corner = point{bg::get<0>(hot_pixels[i]) + expand, bg::get<1>(hot_pixels[i]) + expand};
        std::for_each(
            segments_box_rtree.qbegin(bg::index::intersects(box{min_corner, max_corner})), segments_box_rtree.qend(),
            [&](auto const& val) {
                if (is_point_on_segment(hot_pixels[i], segments[val.second]))
                    segment_pixel_pairs.emplace_back(val.second, i);
            });
    }

    std::vector<edge_t> edges;
    {
        auto [segments_begin, segments_end, pixels] = bucket_sort(
            segment_pixel_pairs, segments.size(), [](auto val) { return val.first; },
            [](auto val) { return val.second; });
        for (std::size_t i = 0; i < segments.size(); i++) {
            auto pixel_begin = std::begin(pixels) + segments_begin[i], pixel_end = std::begin(pixels) + segments_end[i];
            std::sort(pixel_begin, pixel_end, [&](auto pi, auto pj) {
                return less_by_segment{segments[i]}(hot_pixels[pi], hot_pixels[pj]);
            });
            for (; pixel_begin != pixel_end - 1; pixel_begin++) edges.emplace_back(*pixel_begin, *std::next(pixel_begin));
        }
    }
    return {std::move(edges), std::move(hot_pixels)};
}

inline std::vector<edge_with_power_t> edges_to_power(std::vector<edge_t> edges) {
    std::vector<edge_with_power_t> r(edges.size());
    for (std::size_t i = 0; i < r.size(); i++) {
        auto s = edges[i].start, e = edges[i].end;
        if (s < e) r[i] = {s, e, 1};
        else      r[i] = {e, s, -1};
    }
    return r;
}

inline std::vector<edge_with_power_t> unique_edges(std::vector<edge_with_power_t> edges,
                                                     std::size_t num_vertices) {
    auto [begin_loc, end_loc, ordered] = bucket_sort(
        std::move(edges), num_vertices,
        [](const edge_with_power_t& e) { return e.start; },
        [](const edge_with_power_t& e) { return e; });
    std::vector<edge_with_power_t> result;
    for (std::size_t i = 0; i < num_vertices; i++) {
        auto current_begin = std::begin(ordered) + begin_loc[i], current_end = std::begin(ordered) + end_loc[i];
        if (current_begin == current_end) continue;
        std::sort(current_begin, current_end, [](const edge_with_power_t& a, const edge_with_power_t& b) {
            return a.end < b.end;
        });
        for (auto cur = current_begin; cur != current_end; ) {
            int sum = 0;
            auto next = cur;
            while (next != current_end && next->end == cur->end) {
                sum += next->power;
                ++next;
            }
            if (sum != 0)
                result.push_back({cur->start, cur->end, sum});
            cur = next;
        }
    }
    return result;
}

inline auto construct_edges_with_power(auto segments) {
    auto [edges, hot_pixels] = construct_graph(std::move(segments));
    return std::tuple{std::move(hot_pixels),
                      unique_edges(edges_to_power(std::move(edges)), hot_pixels.size())};
}

// ---------------------------------------------------------------------------
// Step 4: chains
// ---------------------------------------------------------------------------

inline chain_build_result build_chains(const std::vector<edge_with_power_t>& sorted_edges,
                                       std::size_t node_num) {
    std::vector<std::size_t> edge_offsets(node_num + 1, 0);
    {
        std::vector<std::size_t> edge_count(node_num, 0);
        for (auto& e : sorted_edges) edge_count[e.start]++;
        std::size_t cur = 0;
        for (std::size_t v = 0; v < node_num; v++) { edge_offsets[v] = cur; cur += edge_count[v]; }
        edge_offsets.back() = cur;
    }

    std::vector<std::uint32_t> out_deg(node_num), in_deg(node_num);
    std::vector<int> out_power(node_num), in_power(node_num);
    for (const auto& e : sorted_edges) {
        out_deg[e.start]++; in_deg[e.end]++; out_power[e.start] = e.power; in_power[e.end] = e.power;
    }

    std::vector<bool> is_end(node_num);
    for (std::size_t i = 0; i < node_num; i++)
        is_end[i] = !(out_deg[i] == 1 && in_deg[i] == 1 && out_power[i] == in_power[i]);

    std::vector<bool> visited(node_num), edge_used(sorted_edges.size());
    std::vector<std::size_t> idx, off{0};
    std::vector<int> powers;
    std::vector<std::size_t> edge_to_chain(sorted_edges.size(), ~0ULL);

    for (std::size_t i = 0; i < node_num; i++) {
        if (!is_end[i]) continue;
        visited[i] = true;
        for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++) {
            if (edge_used[j] || sorted_edges[j].start != i) continue;
            edge_used[j] = true; idx.push_back(i);
            powers.push_back(sorted_edges[j].power); edge_to_chain[j] = off.size() - 1;
            auto cur = sorted_edges[j].end;
            while (!is_end[cur]) {
                visited[cur] = true; idx.push_back(cur);
                std::size_t nj = ~0ULL;
                for (std::size_t k = edge_offsets[cur]; k < edge_offsets[cur + 1]; k++)
                    if (!edge_used[k] && sorted_edges[k].start == cur) { nj = k; break; }
                assert(nj != ~0ULL);
                edge_used[nj] = true; edge_to_chain[nj] = off.size() - 1;
                cur = sorted_edges[nj].end;
            }
            idx.push_back(cur); off.push_back(idx.size());
        }
    }

    for (std::size_t i = 0; i < node_num; i++) {
        if (visited[i]) continue;
        std::size_t sj = ~0ULL;
        for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++)
            if (!edge_used[j] && sorted_edges[j].start == i) { sj = j; break; }
        if (sj == ~0ULL) continue;
        visited[i] = true; edge_used[sj] = true;
        idx.push_back(i); powers.push_back(sorted_edges[sj].power); edge_to_chain[sj] = off.size() - 1;
        auto cur = sorted_edges[sj].end;
        while (i != cur) {
            visited[cur] = true; idx.push_back(cur);
            std::size_t nj = ~0ULL;
            for (std::size_t k = edge_offsets[cur]; k < edge_offsets[cur + 1]; k++)
                if (!edge_used[k] && sorted_edges[k].start == cur) { nj = k; break; }
            assert(nj != ~0ULL);
            edge_used[nj] = true; edge_to_chain[nj] = off.size() - 1;
            cur = sorted_edges[nj].end;
        }
        idx.push_back(cur); off.push_back(idx.size());
    }

    return {std::move(idx), std::move(off), std::move(powers), std::move(edge_to_chain)};
}

// ---------------------------------------------------------------------------
// Step 5-6: angular sort, next/prev, coplanar pairs (sector adjacency)
// ---------------------------------------------------------------------------

inline auto build_half_chain_graph(const chain_build_result& chains,
                           const std::vector<point>& hot_pixels) {
    std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
    std::size_t num_vertices = hot_pixels.size();

    std::vector<half_chain> all;
    all.reserve(num_half_chains);
    for (std::size_t i = 0; i < num_half_chains; i++) all.push_back({i});

    auto [begin_loc, end_loc, sorted_half_chains] = bucket_sort(
        all, num_vertices, [&](half_chain h) { return h.source_node(chains); },
        [](half_chain h) { return h; });

    for (std::size_t v = 0; v < num_vertices; v++) {
        auto vertex_begin = begin_loc[v], vertex_end = end_loc[v];
        if (vertex_end - vertex_begin < 2) continue;
        std::sort(sorted_half_chains.begin() + vertex_begin, sorted_half_chains.begin() + vertex_end,
                  [&](half_chain a, half_chain b) {
                      return less_by_direction(hot_pixels[v],
                                               hot_pixels[a.next_along_source(chains)],
                                               hot_pixels[b.next_along_source(chains)]);
                  });
    }

    std::vector<half_chain> next_half_chain(num_half_chains, {~0ULL});
    std::vector<std::pair<std::size_t, std::size_t>> coplanar;

    for (std::size_t v = 0; v < num_vertices; v++) {
        auto vertex_begin = begin_loc[v], vertex_end = end_loc[v];
        if (vertex_begin == vertex_end) continue;
        assert(vertex_begin + 1 != vertex_end);
        for (auto it = vertex_begin + 1; it < vertex_end; ++it) {
            auto prev = sorted_half_chains[it - 1], cur = sorted_half_chains[it];
            next_half_chain[prev.dual().id] = cur;
            coplanar.emplace_back(cur.id, prev.dual().id);
        }
        auto first = sorted_half_chains[vertex_begin], last = sorted_half_chains[vertex_end - 1];
        next_half_chain[last.dual().id] = first;
        coplanar.emplace_back(first.id, last.dual().id);
    }

    return std::tuple{
        std::move(sorted_half_chains),
        std::vector<std::size_t>(begin_loc.begin(), begin_loc.end()),
        std::vector<std::size_t>(end_loc.begin(), end_loc.end()),
        std::move(next_half_chain), std::move(coplanar)};
}

// ---------------------------------------------------------------------------
// Step 7: ray casting. Returns coplanar pairs from the ray.
// ---------------------------------------------------------------------------

// Cast ray in -x direction from leftmost vertex v.
// Finds the half-chain at v whose right face contains the -x ray, then finds the nearest
// half-chain hit by the ray. These two half-chains are coplanar (share the exterior face).
// If the ray hits nothing, the half-chain at v is the exterior face itself (w=0).
inline std::vector<std::pair<std::size_t, std::size_t>> cast_ray_minus_x(
    std::size_t v,
    const std::vector<point>& hot_pixels,
    const chain_build_result& chains,
    const std::vector<half_chain>& sorted_half_chains,
    const std::vector<std::size_t>& half_chain_begin,
    const std::vector<std::size_t>& half_chain_end) {

    auto range_begin = half_chain_begin[v], range_end = half_chain_end[v];
    if (range_begin == range_end) return {};

    // Step 1: Find the half-chain at v whose right face contains the -x ray direction.
    // The -x ray lies in the sector between prev and cur in CCW order.
    // The sector = right(cur), so cur is the half-chain we want.
    // cur is the first half-chain whose direction is strictly greater than -x in CCW order.
    // If no half-chain direction is greater than -x, the ray wraps around,
    // so cur = the first half-chain (wrapping).
    point vertex_point = hot_pixels[v];
    point ray_point{bg::get<0>(vertex_point) - 1, bg::get<1>(vertex_point)};
    half_chain vertex_half_chain;
    bool found = false;
    for (auto it = range_begin; it < range_end; it++) {
        auto h = sorted_half_chains[it];
        auto half_chain_point = hot_pixels[h.next_along_source(chains)];
        if (less_by_direction(vertex_point, ray_point, half_chain_point)) {
            // ray_dir < half_chain_dir (CCW) — first departure after the ray
            vertex_half_chain = h;
            found = true;
            break;
        }
    }
    if (!found) vertex_half_chain = sorted_half_chains[range_begin];

    // Step 2: Cast ray in -x direction, find the nearest intersected half-chain.
    // Simulation of simplicity: the ray is at y + ε (infinitesimal +dy).
    // This ensures the ray doesn't pass through vertices, so only one half-chain
    // can be hit at any intersection point.
    int64_t ray_y = bg::get<1>(vertex_point);
    int min_x = bg::get<0>(vertex_point);
    int64_t best_x = std::numeric_limits<int64_t>::min();
    std::size_t hit_id = ~0ULL;
    int64_t hit_dy = 0;
    int64_t best_dx = 0;
    int64_t hit_y_low = 0;

    auto try_edge = [&](int64_t ix, int64_t dx, int64_t dy, int64_t y_low, std::size_t half_chain_id) {
        if (ix >= min_x) return;
        bool better = false;
        if (ix > best_x) {
            better = true;
        } else if (ix == best_x) {
            // Compare dx/dy ratio. Edge with larger dx/dy is hit first by
            // the perturbed ray (ix shifts by ε·dx/dy).
            long double slope_a = (long double)dx / dy;
            long double slope_b = (long double)best_dx / hit_dy;
            if (slope_a > slope_b) {
                better = true;
            } else if (slope_a == slope_b) {
                // Equal slopes: use higher-order SoS.
                // Compare parametric position along edge: (ray_y - y_low) / |dy|.
                // Edge hit closer to its bottom endpoint wins.
                int64_t da = ray_y - y_low;         // ≥ 0 since y_low ≤ ray_y
                int64_t db = ray_y - hit_y_low;
                int64_t abs_dy = dy > 0 ? dy : -dy;
                int64_t abs_hdy = hit_dy > 0 ? hit_dy : -hit_dy;
                if (da * abs_hdy < db * abs_dy) {
                    better = true;
                } else if (da * abs_hdy == db * abs_dy && abs_dy > abs_hdy) {
                    // Both hit at same parametric position (e.g. both start at ray_y).
                    // Steeper edge wins (smaller ε/dy for the perturbed ray).
                    better = true;
                }
            }
        }
        if (better) {
            best_x = ix; best_dx = dx; hit_dy = dy; hit_id = half_chain_id;
            hit_y_low = y_low;
        }
    };

    for (auto h : sorted_half_chains) {
        auto& idx = chains.indices;
        auto& off = chains.offsets;
        std::size_t chain_idx = h.chain_id();
        auto chain_begin_idx = off[chain_idx], chain_end_idx = off[chain_idx + 1];

        if (h.is_forward()) {
            for (std::size_t k = chain_begin_idx; k + 1 < chain_end_idx; k++) {
                int64_t y1 = bg::get<1>(hot_pixels[idx[k]]);
                int64_t y2 = bg::get<1>(hot_pixels[idx[k + 1]]);
                if (y1 <= ray_y && ray_y < y2) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k + 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y1, h.id);
                } else if (y2 <= ray_y && ray_y < y1) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k + 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y2, h.id);
                }
            }
        } else {
            for (std::size_t k = chain_end_idx - 1; k > chain_begin_idx; k--) {
                int64_t y1 = bg::get<1>(hot_pixels[idx[k]]);
                int64_t y2 = bg::get<1>(hot_pixels[idx[k - 1]]);
                if (y1 <= ray_y && ray_y < y2) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k - 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y1, h.id);
                } else if (y2 <= ray_y && ray_y < y1) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k - 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y2, h.id);
                }
            }
        }
    }

    // Step 3: Build the coplanar pair (RIGHT-face convention).
    // hit_side is the half-chain whose right face is the approach face of the ray.
    // dy > 0 (edge going up): right face = east, ray approaches from east → hit_id
    // dy < 0 (edge going down): right face = west, ray approaches from east → dual
    // If no hit: vertex_half_chain IS the exterior face.
    std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;
    if (hit_id != ~0ULL) {
        std::size_t hit_side = (hit_dy > 0) ? hit_id : (hit_id ^ 1);
        ray_pairs.emplace_back(vertex_half_chain.id, hit_side);
    } else {
        ray_pairs.emplace_back(vertex_half_chain.id, vertex_half_chain.id);
    }
    return ray_pairs;
}

// Find exterior face info per connected component.
// Returns: the half-chain id in the exterior face (per component),
//   and all ray-coplanar pairs.
inline auto find_exterior(const chain_build_result& chains,
                          const std::vector<point>& hot_pixels,
                          const std::vector<half_chain>& sorted_half_chains,
                          const std::vector<std::size_t>& half_chain_begin,
                          const std::vector<std::size_t>& half_chain_end) {
    std::size_t num_vertices = hot_pixels.size();

    // Vertex connected components
    std::vector<std::pair<std::size_t, std::size_t>> vertex_edges;
    for (auto h : sorted_half_chains)
        vertex_edges.emplace_back(h.source_node(chains), h.target_node(chains));
    for (std::size_t c = 0; c + 1 < chains.offsets.size(); c++) {
        auto chain_begin_idx = chains.offsets[c], chain_end_idx = chains.offsets[c + 1];
        for (std::size_t k = chain_begin_idx; k + 1 < chain_end_idx; k++)
            vertex_edges.emplace_back(chains.indices[k], chains.indices[k + 1]);
    }

    auto component_id = connected_components(num_vertices, vertex_edges);

    // Group vertices by component
    std::size_t num_components = 0;
    for (auto c : component_id) num_components = std::max(num_components, c + 1);
    std::vector<std::vector<std::size_t>> vertex_components(num_components);
    for (std::size_t v = 0; v < num_vertices; v++)
        vertex_components[component_id[v]].push_back(v);

    std::vector<std::size_t> exterior_half_chains;
    std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;

    for (auto& vertex_component : vertex_components) {
        std::size_t leftmost_vertex = ~0ULL;
        int min_x = std::numeric_limits<int>::max();
        for (auto vertex : vertex_component) {
            if (half_chain_begin[vertex] == half_chain_end[vertex]) continue;
            int x = bg::get<0>(hot_pixels[vertex]);
            if (x < min_x) { min_x = x; leftmost_vertex = vertex; }
        }
        if (leftmost_vertex == ~0ULL) continue;

        auto ray_pairs_result = cast_ray_minus_x(leftmost_vertex, hot_pixels, chains, sorted_half_chains, half_chain_begin, half_chain_end);
        for (auto& p : ray_pairs_result) {
            if (p.first == p.second) {
                exterior_half_chains.push_back(p.first);
            }
            ray_pairs.push_back(std::move(p));
        }
    }

    return std::tuple{std::move(exterior_half_chains), std::move(ray_pairs)};
}

// ---------------------------------------------------------------------------
// Step 8: compute winding numbers.
// Coplanarity sources: sector adjacency + ray.
// Propagation: for consecutive a,b at vertex CCW:
//   W(face(b)) = W(face(a)) + power(b)
// Returns: winding for each half-chain.
// ---------------------------------------------------------------------------

inline auto compute_winding(
    const chain_build_result& chains,
    const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
    const std::vector<std::size_t>& exterior_half_chains) {

    std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;

    // Collect edges: coplanar (weight 0), ray (weight 0), dual (weight -power)
    std::vector<edge_with_power_t> edges;
    auto add = [&](std::size_t a, std::size_t b, int w) {
        edges.emplace_back(a, b, w);
        edges.emplace_back(b, a, -w);
    };
    for (auto [a, b] : coplanar_pairs) add(a, b, 0);
    for (auto [a, b] : ray_pairs) add(a, b, 0);
    for (std::size_t i = 0; i < num_half_chains; i++)
        add(i, i ^ 1, -half_chain{i}.power(chains));

    auto [adj_begin, adj_end, adjacency] = bucket_sort(
        edges, num_half_chains,
        [](const edge_with_power_t& e) { return e.start; },
        [](const edge_with_power_t& e) { return std::pair{e.end, e.power}; });

    // DFS propagation from exterior seeds
    constexpr int UNKNOWN = std::numeric_limits<int>::max() / 2;
    std::vector<int> winding(num_half_chains, UNKNOWN);
    for (auto exterior : exterior_half_chains) winding[exterior] = 0;

    std::vector<std::size_t> stack(exterior_half_chains.begin(), exterior_half_chains.end());
    while (!stack.empty()) {
        auto u = stack.back(); stack.pop_back();
        for (auto i = adj_begin[u]; i < adj_end[u]; i++) {
            auto [v, diff] = adjacency[i];
            if (winding[v] == UNKNOWN) {
                winding[v] = winding[u] + diff;
                stack.push_back(v);
            }
        }
    }
    return winding;
}

// ---------------------------------------------------------------------------
// Step 9-10: build polygons face-by-face
// ---------------------------------------------------------------------------

inline multi_polygon build_output(const chain_build_result& chains,
                                   const std::vector<point>& hot_pixels,
                                   std::vector<half_chain> next_half_chain,
                                   const std::vector<int>& winding,
                                   const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
                                   const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
                                   auto filter_fn) {
    std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
    std::size_t num_chains = chains.offsets.size() - 1;

    // 1. Which half-chains survive
    std::vector<bool> survive(num_half_chains);
    for (std::size_t i = 0; i < num_half_chains; i++)
        survive[i] = filter_fn(winding[i]);

    // 2. Dual cancellation: both fwd and rev survive → both dead
    std::vector<bool> dead(num_half_chains, false);
    for (std::size_t chain_idx = 0; chain_idx < num_chains; chain_idx++) {
        std::size_t forward_id = 2 * chain_idx, reverse_id = 2 * chain_idx + 1;
        if (survive[forward_id] && survive[reverse_id])
            dead[forward_id] = dead[reverse_id] = true;
    }

    // 3. Update next_half_chain via indirection: for each surviving half-chain whose next
    //    points to a dead half-chain, follow next_half_chain[dead.dual()] until a live half-chain
    //    is found. This terminates because at least one half-chain per vertex cycle
    //    survives (not all faces can be cancelled).
    for (std::size_t i = 0; i < num_half_chains; i++) {
        if (!survive[i] || dead[i]) continue;
        while (dead[next_half_chain[i].id])
            next_half_chain[i] = next_half_chain[next_half_chain[i].dual().id];
    }

    // 4. Build connected components: coplanar + ray + dual cancellation → same face
    std::vector<std::pair<std::size_t, std::size_t>> face_edges;
    face_edges.insert(face_edges.end(), coplanar_pairs.begin(), coplanar_pairs.end());
    face_edges.insert(face_edges.end(), ray_pairs.begin(), ray_pairs.end());
    for (std::size_t chain_idx = 0; chain_idx < num_chains; chain_idx++) {
        std::size_t forward_id = 2 * chain_idx, reverse_id = 2 * chain_idx + 1;
        if (dead[forward_id]) face_edges.emplace_back(forward_id, reverse_id);
    }
    auto component_id = connected_components(num_half_chains, face_edges);

    // 5. Group surviving non-dead half-chains by face via bucket_sort
    std::vector<std::pair<std::size_t, std::size_t>> half_chain_to_face;
    for (std::size_t i = 0; i < num_half_chains; i++)
        if (!dead[i] && survive[i])
            half_chain_to_face.emplace_back(component_id[i], i);

    std::size_t num_faces = 0;
    for (auto c : component_id) num_faces = std::max(num_faces, c + 1);
    auto [face_begin, face_end, face_data] = bucket_sort(
        half_chain_to_face, num_faces,
        [](const std::pair<std::size_t, std::size_t>& p) { return p.first; },
        [](const std::pair<std::size_t, std::size_t>& p) { return p.second; });

    // 6. For each face, trace its half-chains into a polygon
    multi_polygon result;

    for (std::size_t f = 0; f < num_faces; f++) {
        auto face_data_begin = face_begin[f], face_data_end = face_end[f];
        if (face_data_begin == face_data_end) continue;
        std::size_t half_chain_count = face_data_end - face_data_begin;

        boost::container::flat_map<std::size_t, std::size_t> position_map;
        for (std::size_t j = 0; j < half_chain_count; j++)
            position_map[face_data[face_data_begin + j]] = j;

        std::vector<bool> done(half_chain_count, false);
        std::vector<ring> rings;

        for (std::size_t j = 0; j < half_chain_count; j++) {
            if (done[j]) continue;

            ring current_ring;
            std::size_t current_id = face_data[face_data_begin + j];
            do {
                auto h = half_chain{current_id};
                auto chain_idx = h.chain_id();
                auto chain_begin = chains.offsets[chain_idx], chain_end = chains.offsets[chain_idx + 1];
                if (h.is_forward()) {
                    for (std::size_t k = chain_begin; k < chain_end - 1; k++)
                        current_ring.push_back(hot_pixels[chains.indices[k]]);
                } else {
                    for (std::size_t k = chain_end - 1; k > chain_begin; k--)
                        current_ring.push_back(hot_pixels[chains.indices[k]]);
                }
                auto it = position_map.find(current_id);
                if (it != position_map.end()) done[it->second] = true;
                current_id = next_half_chain[current_id].id;
            } while (current_id != face_data[face_data_begin + j]);

            if (!current_ring.empty()) current_ring.push_back(current_ring.front());
            rings.push_back(std::move(current_ring));
        }

        if (rings.empty()) continue;

        std::size_t outer_index = 0;
        {
            int min_x = std::numeric_limits<int>::max();
            for (std::size_t ring_index = 0; ring_index < rings.size(); ring_index++) {
                for (const auto& point_value : rings[ring_index]) {
                    int x = bg::get<0>(point_value);
                    if (x < min_x) { min_x = x; outer_index = ring_index; }
                }
            }
        }

        polygon polygon_result;
        polygon_result.outer() = std::move(rings[outer_index]);
        for (std::size_t ring_index = 0; ring_index < rings.size(); ring_index++) {
            if (ring_index == outer_index || rings[ring_index].empty()) continue;
            polygon_result.inners().push_back(std::move(rings[ring_index]));
        }

        result.push_back(std::move(polygon_result));
    }

    return result;
}

// ---------------------------------------------------------------------------
// Main pipeline
// ---------------------------------------------------------------------------

inline auto run_pipeline(auto segments, auto filter) {
    // Phase 1: edges → chains → half-chain graph
    auto [hot_pixels, sorted_edges] = construct_edges_with_power(std::move(segments));

    auto chains = build_chains(sorted_edges, hot_pixels.size());

    auto [sorted_half_chains, half_chain_begin, half_chain_end, next_half_chain, coplanar] =
        build_half_chain_graph(chains, hot_pixels);

    // Phase 2: ray casting → ray coplanar pairs
    auto [exterior_half_chains, ray_pairs] =
        find_exterior(chains, hot_pixels, sorted_half_chains, half_chain_begin, half_chain_end);

    // Phase 3: winding numbers
    auto winding = compute_winding(chains, coplanar, ray_pairs, exterior_half_chains);

    // Phase 4: build polygons face by face
    return build_output(chains, hot_pixels, std::move(next_half_chain),
                        winding, coplanar, ray_pairs, filter);
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

inline auto collect_segments(auto... polygons) {
    std::vector<segment> segments;
    ((bg::for_each_segment(polygons, [&](const auto& seg) {
        segment s;
        bg::set<0, 0>(s, bg::get<0, 0>(seg)); bg::set<0, 1>(s, bg::get<0, 1>(seg));
        bg::set<1, 0>(s, bg::get<1, 0>(seg)); bg::set<1, 1>(s, bg::get<1, 1>(seg));
        segments.emplace_back(s);
    })), ...);
    return segments;
}

inline auto add(const auto& polygons1, const auto& polygons2) {
    return run_pipeline(collect_segments(polygons1, polygons2), [](int w) { return w > 0; });
}

inline auto intersection(const auto& polygons1, const auto& polygons2) {
    return run_pipeline(collect_segments(polygons1, polygons2), [](int w) { return w > 1; });
}

inline auto self_or(auto r) {
    return run_pipeline(collect_segments(r), [](int w) { return w > 0; });
}

inline auto xor_(const auto& polygons1, const auto& polygons2) {
    return run_pipeline(collect_segments(polygons1, polygons2), [](int w) { return w == 1; });
}

inline auto difference(const auto& polygons1, const auto& polygons2) {
    auto segments = collect_segments(polygons1);
    bg::for_each_segment(polygons2, [&](const auto& seg) {
        segment s;
        bg::set<0, 0>(s, bg::get<1, 0>(seg)); bg::set<0, 1>(s, bg::get<1, 1>(seg));
        bg::set<1, 0>(s, bg::get<0, 0>(seg)); bg::set<1, 1>(s, bg::get<0, 1>(seg));
        segments.emplace_back(s);
    });
    return run_pipeline(std::move(segments), [](int w) { return w > 0; });
}
