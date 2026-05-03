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

inline std::tuple<std::vector<edge_t>, std::vector<point>> construct_graph(auto segs) {
    std::vector<std::pair<box, std::size_t>> boxes(segs.size());
    for (std::size_t i = 0; i < segs.size(); i++)
        boxes[i] = {bg::return_envelope<box>(segs[i]), i};
    bg::index::rtree<std::pair<box, std::size_t>, bg::index::quadratic<128>> segs_box_rtree(
        std::move(boxes));

    std::vector<point> hot_pixels;
    hot_pixels.reserve(segs.size() * 2);
    for (const auto& seg : segs) {
        hot_pixels.emplace_back(bg::get<0, 0>(seg), bg::get<0, 1>(seg));
        hot_pixels.emplace_back(bg::get<1, 0>(seg), bg::get<1, 1>(seg));
    }

    for (std::size_t i = 0; i < segs.size(); i++)
        std::for_each(
            segs_box_rtree.qbegin(bg::index::intersects(boxes[i].first)), segs_box_rtree.qend(),
            [&](auto const& o) {
                if (auto p = get_intersection(segs[i], segs[o.second]))
                    hot_pixels.push_back(p.value());
            });

    {
        std::sort(std::begin(hot_pixels), std::end(hot_pixels), [](auto p1, auto p2) {
            return std::pair{bg::get<0>(p1), bg::get<1>(p1)} <
                   std::pair{bg::get<0>(p2), bg::get<1>(p2)};
        });
        auto last = std::unique(std::begin(hot_pixels), std::end(hot_pixels),
                                [](auto p1, auto p2) { return bg::equals(p1, p2); });
        hot_pixels.erase(last, hot_pixels.end());
        bg::index::rtree<point, bg::index::quadratic<128>> hp_rtree{hot_pixels};
        // Keep hot_pixels in xy order — edges_to_power relies on indices being xy-ordered
        // so that start > end means "against natural direction" consistently.
    }

    std::vector<std::pair<std::size_t, std::size_t>> seg_pixel_pairs;
    for (std::size_t i = 0; i < hot_pixels.size(); i++) {
        constexpr auto expand = 1;
        auto mc = point{bg::get<0>(hot_pixels[i]) - expand, bg::get<1>(hot_pixels[i]) - expand};
        auto Mc = point{bg::get<0>(hot_pixels[i]) + expand, bg::get<1>(hot_pixels[i]) + expand};
        std::for_each(
            segs_box_rtree.qbegin(bg::index::intersects(box{mc, Mc})), segs_box_rtree.qend(),
            [&](auto const& val) {
                if (is_point_on_segment(hot_pixels[i], segs[val.second]))
                    seg_pixel_pairs.emplace_back(val.second, i);
            });
    }

    std::vector<edge_t> edges;
    {
        auto [segs_begin, segs_end, pixels] = bucket_sort(
            seg_pixel_pairs, segs.size(), [](auto val) { return val.first; },
            [](auto val) { return val.second; });
        for (std::size_t i = 0; i < segs.size(); i++) {
            auto cb = std::begin(pixels) + segs_begin[i], ce = std::begin(pixels) + segs_end[i];
            std::sort(cb, ce, [&](auto pi, auto pj) {
                return less_by_segment{segs[i]}(hot_pixels[pi], hot_pixels[pj]);
            });
            for (; cb != ce - 1; cb++) edges.emplace_back(*cb, *std::next(cb));
        }
    }
    return {std::move(edges), std::move(hot_pixels)};
}

inline std::vector<edge_with_power_t> edges_to_power(std::vector<edge_t> edges) {
    std::vector<edge_with_power_t> r(edges.size());
    for (std::size_t i = 0; i < r.size(); i++) {
        r[i] = {edges[i].start, edges[i].end, 1};
        if (r[i].start > r[i].end) { std::swap(r[i].start, r[i].end); r[i].power = -1; }
    }
    return r;
}

inline std::vector<edge_with_power_t> sort_edges_by_start(std::vector<edge_with_power_t> edges,
                                                           std::size_t nv) {
    auto [bl, el, ordered] = bucket_sort(
        std::move(edges), nv, [](edge_with_power_t e) { return e.start; },
        [](edge_with_power_t e) { return e; });
    for (std::size_t i = 0; i < nv; i++) {
        auto cb = std::begin(ordered) + bl[i], ce = std::begin(ordered) + el[i];
        std::sort(cb, ce, [](auto a, auto b) { return a.end < b.end; });
    }
    return ordered;
}

inline std::vector<edge_with_power_t> unique_edges(std::vector<edge_with_power_t> edges) {
    constexpr auto eq = [](edge_with_power_t a, edge_with_power_t b) {
        return a.start == b.start && a.end == b.end;
    };
    constexpr auto mg = [](edge_with_power_t a, edge_with_power_t b) {
        return edge_with_power_t{a.start, a.end, a.power + b.power};
    };
    auto cur = std::begin(edges), out = cur;
    while (cur != std::end(edges)) {
        auto nxt = not_adjacent_find(cur, std::end(edges), eq);
        *out++ = std::reduce(std::next(cur), nxt, *cur, mg);
        cur = nxt;
    }
    edges.erase(out, std::end(edges));
    return edges;
}

inline auto construct_edges_with_power(auto segs) {
    auto [edges, hp] = construct_graph(std::move(segs));
    auto ewp = unique_edges(sort_edges_by_start(edges_to_power(std::move(edges)), hp.size()));
    std::erase_if(ewp, [](auto e) { return e.power == 0; });
    return std::tuple{std::move(hp), std::move(ewp)};
}

// ---------------------------------------------------------------------------
// Step 4: chains
// ---------------------------------------------------------------------------

inline std::size_t edge_dest(const edge_with_power_t& e) { return e.power > 0 ? e.end : e.start; }
inline std::size_t edge_ts(const edge_with_power_t& e) { return e.power > 0 ? e.start : e.end; }

inline chain_build_result build_chains(const std::vector<edge_with_power_t>& sorted_edges,
                                       std::size_t node_num) {
    std::vector<std::size_t> edge_offsets(node_num + 1, 0);
    {
        std::vector<std::size_t> cnt(node_num, 0);
        for (auto& e : sorted_edges) cnt[edge_ts(e)]++;
        std::size_t cur = 0;
        for (std::size_t v = 0; v < node_num; v++) { edge_offsets[v] = cur; cur += cnt[v]; }
        edge_offsets.back() = cur;
    }

    std::vector<std::uint32_t> out_deg(node_num), in_deg(node_num);
    std::vector<int> out_pow(node_num), in_pow(node_num);
    for (const auto& e : sorted_edges) {
        auto s = edge_ts(e), d = edge_dest(e);
        out_deg[s]++; in_deg[d]++; out_pow[s] = e.power; in_pow[d] = e.power;
    }

    std::vector<bool> is_end(node_num);
    for (std::size_t i = 0; i < node_num; i++)
        is_end[i] = !(out_deg[i] == 1 && in_deg[i] == 1 && out_pow[i] == in_pow[i]);

    std::vector<bool> visited(node_num), edge_used(sorted_edges.size());
    std::vector<std::size_t> idx, off{0};
    std::vector<int> pows;
    std::vector<std::size_t> e2c(sorted_edges.size(), ~0ULL);

    for (std::size_t i = 0; i < node_num; i++) {
        if (!is_end[i]) continue;
        visited[i] = true;
        for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++) {
            if (edge_used[j] || edge_ts(sorted_edges[j]) != i) continue;
            edge_used[j] = true; idx.push_back(i);
            pows.push_back(std::abs(sorted_edges[j].power)); e2c[j] = off.size() - 1;
            auto cur = edge_dest(sorted_edges[j]);
            while (!is_end[cur]) {
                visited[cur] = true; idx.push_back(cur);
                std::size_t nj = ~0ULL;
                for (std::size_t k = edge_offsets[cur]; k < edge_offsets[cur + 1]; k++)
                    if (!edge_used[k] && edge_ts(sorted_edges[k]) == cur) { nj = k; break; }
                assert(nj != ~0ULL);
                edge_used[nj] = true; e2c[nj] = off.size() - 1;
                cur = edge_dest(sorted_edges[nj]);
            }
            idx.push_back(cur); off.push_back(idx.size());
        }
    }

    for (std::size_t i = 0; i < node_num; i++) {
        if (visited[i]) continue;
        std::size_t sj = ~0ULL;
        for (std::size_t j = edge_offsets[i]; j < edge_offsets[i + 1]; j++)
            if (!edge_used[j] && edge_ts(sorted_edges[j]) == i) { sj = j; break; }
        if (sj == ~0ULL) continue;
        visited[i] = true; edge_used[sj] = true;
        idx.push_back(i); pows.push_back(std::abs(sorted_edges[sj].power)); e2c[sj] = off.size() - 1;
        auto cur = edge_dest(sorted_edges[sj]);
        while (i != cur) {
            visited[cur] = true; idx.push_back(cur);
            std::size_t nj = ~0ULL;
            for (std::size_t k = edge_offsets[cur]; k < edge_offsets[cur + 1]; k++)
                if (!edge_used[k] && edge_ts(sorted_edges[k]) == cur) { nj = k; break; }
            assert(nj != ~0ULL);
            edge_used[nj] = true; e2c[nj] = off.size() - 1;
            cur = edge_dest(sorted_edges[nj]);
        }
        idx.push_back(cur); off.push_back(idx.size());
    }

    return {std::move(idx), std::move(off), std::move(pows), std::move(e2c)};
}

// ---------------------------------------------------------------------------
// DSU helper
// ---------------------------------------------------------------------------

struct dsu_t {
    std::vector<std::size_t> p;
    std::size_t find(std::size_t x) {
        std::size_t r = x;
        while (p[r] != r) r = p[r];
        while (x != r) { auto n = p[x]; p[x] = r; x = n; }
        return r;
    }
    void unite(std::size_t a, std::size_t b) {
        a = find(a); b = find(b);
        if (a != b) p[a] = b;
    }
    void compress() {
        for (std::size_t i = 0; i < p.size(); i++) p[i] = find(i);
    }
};

inline dsu_t make_dsu(std::size_t n) {
    dsu_t d; d.p.resize(n);
    for (std::size_t i = 0; i < n; i++) d.p[i] = i;
    return d;
}

// ---------------------------------------------------------------------------
// Step 5-6: angular sort, next/prev, coplanar pairs (sector adjacency)
// ---------------------------------------------------------------------------

inline auto build_hc_graph(const chain_build_result& chains,
                           const std::vector<point>& hot_pixels) {
    std::size_t num_hcs = (chains.offsets.size() - 1) * 2;
    std::size_t nv = hot_pixels.size();

    std::vector<half_chain> all;
    all.reserve(num_hcs);
    for (std::size_t i = 0; i < num_hcs; i++) all.push_back({i});

    auto [begin_loc, end_loc, sorted_hcs] = bucket_sort(
        all, nv, [&](half_chain hc) { return hc.source_node(chains); },
        [](half_chain hc) { return hc; });

    for (std::size_t v = 0; v < nv; v++) {
        auto cb = begin_loc[v], ce = end_loc[v];
        if (ce - cb < 2) continue;
        std::sort(sorted_hcs.begin() + cb, sorted_hcs.begin() + ce,
                  [&](half_chain a, half_chain b) {
                      return less_by_direction(hot_pixels[v],
                                               hot_pixels[a.next_along_source(chains)],
                                               hot_pixels[b.next_along_source(chains)]);
                  });
    }

    std::vector<half_chain> next_hc(num_hcs, {~0ULL});
    std::vector<std::pair<std::size_t, std::size_t>> coplanar;

    for (std::size_t v = 0; v < nv; v++) {
        auto cb = begin_loc[v], ce = end_loc[v];
        if (cb == ce) continue;
        assert(cb + 1 != ce);
        for (auto it = cb + 1; it < ce; ++it) {
            auto prev = sorted_hcs[it - 1], cur = sorted_hcs[it];
            next_hc[prev.dual().id] = cur;
            coplanar.emplace_back(cur.id, prev.dual().id);
        }
        auto first = sorted_hcs[cb], last = sorted_hcs[ce - 1];
        next_hc[last.dual().id] = first;
        coplanar.emplace_back(first.id, last.dual().id);
    }

    return std::tuple{
        std::move(sorted_hcs),
        std::vector<std::size_t>(begin_loc.begin(), begin_loc.end()),
        std::vector<std::size_t>(end_loc.begin(), end_loc.end()),
        std::move(next_hc), std::move(coplanar)};
}

// ---------------------------------------------------------------------------
// Step 7: ray casting. Returns coplanar pairs from the ray.
// ---------------------------------------------------------------------------

// Cast ray in -x direction from leftmost vertex v.
// Finds the HC at v whose right face contains the -x ray, then finds the nearest
// HC hit by the ray. These two HCs are coplanar (share the exterior face).
// If the ray hits nothing, the HC at v is the exterior face itself (w=0).
inline std::vector<std::pair<std::size_t, std::size_t>> cast_ray_minus_x(
    std::size_t v,
    const std::vector<point>& hot_pixels,
    const chain_build_result& chains,
    const std::vector<half_chain>& sorted_hcs,
    const std::vector<std::size_t>& hc_begin,
    const std::vector<std::size_t>& hc_end) {

    auto beg = hc_begin[v], end = hc_end[v];
    if (beg == end) return {};

    // Step 1: Find the HC at v whose right face contains the -x ray direction.
    // The -x ray lies in the sector between prev and cur in CCW order.
    // The sector = right(cur), so cur is the HC we want.
    // cur is the first HC whose direction is strictly greater than -x in CCW order.
    // If no HC direction is greater than -x, the ray wraps around,
    // so cur = the first HC (wrapping).
    point v_pt = hot_pixels[v];
    point ray_pt{bg::get<0>(v_pt) - 1, bg::get<1>(v_pt)};
    half_chain hc_vertex;
    bool found = false;
    for (auto it = beg; it < end; it++) {
        auto hc = sorted_hcs[it];
        auto hc_pt = hot_pixels[hc.next_along_source(chains)];
        if (less_by_direction(v_pt, ray_pt, hc_pt)) {
            // ray_dir < hc_dir (CCW) — first departure after the ray
            hc_vertex = hc;
            found = true;
            break;
        }
    }
    if (!found) hc_vertex = sorted_hcs[beg];

    // Step 2: Cast ray in -x direction, find the nearest intersected HC.
    // Simulation of simplicity: the ray is at y + ε (infinitesimal +dy).
    // This ensures the ray doesn't pass through vertices, so only one HC
    // can be hit at any intersection point.
    int64_t ray_y = bg::get<1>(v_pt);
    int min_x = bg::get<0>(v_pt);
    int64_t best_x = std::numeric_limits<int64_t>::min();
    std::size_t hit_id = ~0ULL;
    int64_t hit_dy = 0;
    int64_t best_dx = 0;
    int64_t hit_y_low = 0;

    auto try_edge = [&](int64_t ix, int64_t dx, int64_t dy, int64_t y_low, std::size_t dc_id) {
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
            best_x = ix; best_dx = dx; hit_dy = dy; hit_id = dc_id;
            hit_y_low = y_low;
        }
    };

    for (auto hc : sorted_hcs) {
        auto& idx = chains.indices;
        auto& off = chains.offsets;
        std::size_t cid = hc.chain_id();
        auto cb = off[cid], ce = off[cid + 1];

        if (hc.is_forward()) {
            for (std::size_t k = cb; k + 1 < ce; k++) {
                int64_t y1 = bg::get<1>(hot_pixels[idx[k]]);
                int64_t y2 = bg::get<1>(hot_pixels[idx[k + 1]]);
                if (y1 <= ray_y && ray_y < y2) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k + 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y1, hc.id);
                } else if (y2 <= ray_y && ray_y < y1) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k + 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y2, hc.id);
                }
            }
        } else {
            for (std::size_t k = ce - 1; k > cb; k--) {
                int64_t y1 = bg::get<1>(hot_pixels[idx[k]]);
                int64_t y2 = bg::get<1>(hot_pixels[idx[k - 1]]);
                if (y1 <= ray_y && ray_y < y2) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k - 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y1, hc.id);
                } else if (y2 <= ray_y && ray_y < y1) {
                    int64_t x1 = bg::get<0>(hot_pixels[idx[k]]);
                    int64_t x2 = bg::get<0>(hot_pixels[idx[k - 1]]);
                    long double t = (long double)(ray_y - y1) / (y2 - y1);
                    int64_t ix = (int64_t)(x1 + t * (x2 - x1));
                    try_edge(ix, x2 - x1, y2 - y1, y2, hc.id);
                }
            }
        }
    }

    // Step 3: Build the coplanar pair (RIGHT-face convention).
    // hit_side is the HC whose right face is the approach face of the ray.
    // dy > 0 (edge going up): right face = east, ray approaches from east → hit_id
    // dy < 0 (edge going down): right face = west, ray approaches from east → dual
    // If no hit: hc_vertex IS the exterior face.
    std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;
    if (hit_id != ~0ULL) {
        std::size_t hit_side = (hit_dy > 0) ? hit_id : (hit_id ^ 1);
        ray_pairs.emplace_back(hc_vertex.id, hit_side);
    } else {
        ray_pairs.emplace_back(hc_vertex.id, hc_vertex.id);
    }
    return ray_pairs;
}

// Find exterior face info per connected component.
// Returns: the HC id in the exterior face (per component),
//   and all ray-coplanar pairs.
inline auto find_exterior(const chain_build_result& chains,
                          const std::vector<point>& hot_pixels,
                          const std::vector<half_chain>& sorted_hcs,
                          const std::vector<std::size_t>& hc_begin,
                          const std::vector<std::size_t>& hc_end) {
    std::size_t nv = hot_pixels.size();

    // Vertex components: union HC endpoints AND consecutive chain vertices
    auto dsu_v = make_dsu(nv);
    for (auto hc : sorted_hcs)
        dsu_v.unite(hc.source_node(chains), hc.target_node(chains));
    for (std::size_t c = 0; c + 1 < chains.offsets.size(); c++) {
        auto cb = chains.offsets[c], ce = chains.offsets[c + 1];
        for (std::size_t k = cb; k + 1 < ce; k++)
            dsu_v.unite(chains.indices[k], chains.indices[k + 1]);
    }
    dsu_v.compress();

    std::vector<std::vector<std::size_t>> comps(nv);
    for (std::size_t v = 0; v < nv; v++) comps[dsu_v.p[v]].push_back(v);

    std::vector<std::size_t> ext_hcs;
    std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;

    for (std::size_t c = 0; c < nv; c++) {
        if (comps[c].empty()) continue;

        std::size_t lmv = ~0ULL;
        int min_x = std::numeric_limits<int>::max();
        for (auto vv : comps[c]) {
            if (hc_begin[vv] == hc_end[vv]) continue;
            int x = bg::get<0>(hot_pixels[vv]);
            if (x < min_x) { min_x = x; lmv = vv; }
        }
        if (lmv == ~0ULL) continue;

        auto rp = cast_ray_minus_x(lmv, hot_pixels, chains, sorted_hcs, hc_begin, hc_end);
        for (auto& p : rp) {
            if (p.first == p.second) {
                ext_hcs.push_back(p.first);
            }
            ray_pairs.push_back(std::move(p));
        }
    }

    return std::tuple{std::move(ext_hcs), std::move(ray_pairs)};
}

// ---------------------------------------------------------------------------
// Step 8: compute winding numbers.
// Coplanarity sources: sector adjacency + ray.
// Propagation: for consecutive a,b at vertex CCW:
//   W(face(b)) = W(face(a)) + power(b)
// Returns: hc_winding for each HC.
// ---------------------------------------------------------------------------

inline auto compute_hc_winding(
    const chain_build_result& chains,
    const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
    const std::vector<std::size_t>& exterior_hcs) {

    std::size_t num_hcs = (chains.offsets.size() - 1) * 2;

    // Collect edges: coplanar (weight 0), ray (weight 0), dual (weight -power)
    std::vector<std::tuple<std::size_t, std::size_t, int>> edges;
    auto add = [&](std::size_t a, std::size_t b, int w) {
        edges.emplace_back(a, b, w);
        edges.emplace_back(b, a, -w);
    };
    for (auto [a, b] : coplanar_pairs) add(a, b, 0);
    for (auto [a, b] : ray_pairs) add(a, b, 0);
    for (std::size_t i = 0; i < num_hcs; i++)
        add(i, i ^ 1, -half_chain{i}.power(chains));

    auto [adj_begin, adj_end, adj] = bucket_sort(
        edges, num_hcs,
        [](auto& e) { return std::get<0>(e); },
        [](auto& e) { return std::pair{std::get<1>(e), std::get<2>(e)}; });

    // DFS propagation from exterior seeds
    constexpr int UNK = std::numeric_limits<int>::max() / 2;
    std::vector<int> winding(num_hcs, UNK);
    for (auto ext : exterior_hcs) winding[ext] = 0;

    std::vector<std::size_t> stack(exterior_hcs.begin(), exterior_hcs.end());
    while (!stack.empty()) {
        auto u = stack.back(); stack.pop_back();
        for (auto i = adj_begin[u]; i < adj_end[u]; i++) {
            auto [v, diff] = adj[i];
            if (winding[v] == UNK) {
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
                                   std::vector<half_chain> next_hc,
                                   const std::vector<int>& hc_winding,
                                   const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
                                   const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
                                   auto filter_fn) {
    std::size_t num_hcs = (chains.offsets.size() - 1) * 2;
    std::size_t num_chains = chains.offsets.size() - 1;

    // 1. Which HCs survive
    std::vector<bool> survive(num_hcs);
    for (std::size_t i = 0; i < num_hcs; i++)
        survive[i] = filter_fn(hc_winding[i]);

    // 2. Dual cancellation: both fwd and rev survive → both dead
    std::vector<bool> dead(num_hcs, false);
    for (std::size_t cid = 0; cid < num_chains; cid++) {
        std::size_t fwd = 2 * cid, rev = 2 * cid + 1;
        if (survive[fwd] && survive[rev])
            dead[fwd] = dead[rev] = true;
    }

    // 3. Update next_hc via indirection: for each surviving HC whose next
    //    points to a dead HC, follow next_hc[dead.dual()] until a live HC
    //    is found. This terminates because at least one HC per vertex cycle
    //    survives (not all faces can be cancelled).
    for (std::size_t i = 0; i < num_hcs; i++) {
        if (!survive[i] || dead[i]) continue;
        while (dead[next_hc[i].id])
            next_hc[i] = next_hc[next_hc[i].dual().id];
    }

    // 4. Build DSU: coplanar pairs + ray pairs + dual cancellation merges
    auto dsu = make_dsu(num_hcs);
    for (auto [a, b] : coplanar_pairs) dsu.unite(a, b);
    for (auto [a, b] : ray_pairs) dsu.unite(a, b);
    for (std::size_t cid = 0; cid < num_chains; cid++) {
        std::size_t fwd = 2 * cid, rev = 2 * cid + 1;
        if (dead[fwd]) dsu.unite(fwd, rev);
    }
    dsu.compress();

    // 5. Group surviving non-dead HCs by face
    std::vector<std::vector<std::size_t>> face_hcs(num_hcs);
    for (std::size_t i = 0; i < num_hcs; i++) {
        if (!dead[i] && survive[i])
            face_hcs[dsu.p[i]].push_back(i);
    }

    // 4. For each face, trace its HCs into a polygon
    multi_polygon result;

    for (std::size_t f = 0; f < num_hcs; f++) {
        auto& dcs = face_hcs[f];
        if (dcs.empty()) continue;

        // Map dc_id → position in dcs
        boost::container::flat_map<std::size_t, std::size_t> pos;
        for (std::size_t j = 0; j < dcs.size(); j++) pos[dcs[j]] = j;

        std::vector<bool> done(dcs.size(), false);
        std::vector<ring> rings;

        for (std::size_t j = 0; j < dcs.size(); j++) {
            if (done[j]) continue;

            ring r;
            std::size_t cur_id = dcs[j];
            do {
                auto hc = half_chain{cur_id};
                auto cid = hc.chain_id();
                auto cb = chains.offsets[cid], ce = chains.offsets[cid + 1];
                if (hc.is_forward()) {
                    for (std::size_t k = cb; k < ce - 1; k++)
                        r.push_back(hot_pixels[chains.indices[k]]);
                } else {
                    for (std::size_t k = ce - 1; k > cb; k--)
                        r.push_back(hot_pixels[chains.indices[k]]);
                }
                auto it = pos.find(cur_id);
                if (it != pos.end()) done[it->second] = true;
                cur_id = next_hc[cur_id].id;
            } while (cur_id != dcs[j]);

            // Close ring (boost::geometry requires closed rings)
            if (!r.empty()) r.push_back(r.front());

            rings.push_back(std::move(r));
        }

        if (rings.empty()) continue;

        // Find outer ring: the one containing the leftmost point
        std::size_t outer_i = 0;
        {
            int min_x = std::numeric_limits<int>::max();
            for (std::size_t r = 0; r < rings.size(); r++) {
                for (const auto& pt : rings[r]) {
                    int x = bg::get<0>(pt);
                    if (x < min_x) { min_x = x; outer_i = r; }
                }
            }
        }

        polygon poly;
        poly.outer() = std::move(rings[outer_i]);
        for (std::size_t r = 0; r < rings.size(); r++) {
            if (r == outer_i || rings[r].empty()) continue;
            poly.inners().push_back(std::move(rings[r]));
        }

        result.push_back(std::move(poly));
    }

    return result;
}

// ---------------------------------------------------------------------------
// Main pipeline
// ---------------------------------------------------------------------------

inline auto run_pipeline(auto segs, auto filter) {
    // Phase 1: edges → chains → HC graph
    auto [hot_pixels, ewp] = construct_edges_with_power(std::move(segs));

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

    auto chains = build_chains(sorted_edges, hot_pixels.size());

    auto [sorted_hcs, hc_begin, hc_end, next_hc, coplanar] =
        build_hc_graph(chains, hot_pixels);

    // Phase 2: ray casting → ray coplanar pairs
    auto [ext_hcs, ray_pairs] =
        find_exterior(chains, hot_pixels, sorted_hcs, hc_begin, hc_end);

    // Phase 3: winding numbers
    auto hc_winding = compute_hc_winding(chains, coplanar, ray_pairs, ext_hcs);

    // Phase 4: build polygons face by face
    return build_output(chains, hot_pixels, std::move(next_hc),
                        hc_winding, coplanar, ray_pairs, filter);
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

inline auto collect_segments(auto... ps) {
    std::vector<segment> segs;
    ((bg::for_each_segment(ps, [&](const auto& seg) {
        segment s;
        bg::set<0, 0>(s, bg::get<0, 0>(seg)); bg::set<0, 1>(s, bg::get<0, 1>(seg));
        bg::set<1, 0>(s, bg::get<1, 0>(seg)); bg::set<1, 1>(s, bg::get<1, 1>(seg));
        segs.emplace_back(s);
    })), ...);
    return segs;
}

inline auto add(const auto& ps1, const auto& ps2) {
    return run_pipeline(collect_segments(ps1, ps2), [](int w) { return w > 0; });
}

inline auto intersection(const auto& ps1, const auto& ps2) {
    return run_pipeline(collect_segments(ps1, ps2), [](int w) { return w > 1; });
}

inline auto self_or(auto r) {
    return run_pipeline(collect_segments(r), [](int w) { return w > 0; });
}

inline auto xor_(const auto& ps1, const auto& ps2) {
    return run_pipeline(collect_segments(ps1, ps2), [](int w) { return w == 1; });
}

inline auto difference(const auto& ps1, const auto& ps2) {
    auto segs = collect_segments(ps1);
    bg::for_each_segment(ps2, [&](const auto& seg) {
        segment s;
        bg::set<0, 0>(s, bg::get<1, 0>(seg)); bg::set<0, 1>(s, bg::get<1, 1>(seg));
        bg::set<1, 0>(s, bg::get<0, 0>(seg)); bg::set<1, 1>(s, bg::get<0, 1>(seg));
        segs.emplace_back(s);
    });
    return run_pipeline(std::move(segs), [](int w) { return w > 0; });
}
