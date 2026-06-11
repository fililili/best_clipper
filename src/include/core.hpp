#pragma once

#include <algorithm>
#include <boost/geometry.hpp>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <vector>

#ifdef _MSC_VER
#include <boost/multiprecision/cpp_int.hpp>
#endif

namespace best_clipper {

#ifdef _MSC_VER
using int128_t = boost::multiprecision::int128_t;
#else
using int128_t = __int128;
#endif

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<int32_t>;
using segment = bg::model::segment<point>;
using box = bg::model::box<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

inline auto bucket_sort(auto vec, auto bucket_size, auto get_bucket, auto get_left) {
    std::vector<std::size_t> times(bucket_size);
    for (auto val : vec) times[get_bucket(val)]++;
    std::vector<std::size_t> begin_location(times.size());
    std::exclusive_scan(std::begin(times), std::end(times), std::begin(begin_location), std::size_t{0});
    std::vector<std::invoke_result_t<decltype(get_left), typename decltype(vec)::value_type>> left(
        vec.size());
    auto current_location{begin_location};
    for (auto val : vec) left[current_location[get_bucket(val)]++] = get_left(val);
    return std::tuple{std::move(begin_location), std::move(current_location), std::move(left)};
}

// Undirected graph connected components via DFS. O(n + m).
std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>>& edges);

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

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
// Edge construction (intersections, split, dedup, assign power)
// ---------------------------------------------------------------------------

std::tuple<std::vector<point>, std::vector<edge_with_power_t>>
construct_edges_with_power(const std::vector<point>& points, const std::vector<std::size_t>& offsets);

// ---------------------------------------------------------------------------
// Chains
// ---------------------------------------------------------------------------

chain_build_result build_chains(const std::vector<edge_with_power_t>& sorted_edges,
                                 std::size_t node_num);

// ---------------------------------------------------------------------------
// Half-chain graph (angular sort, next pointers, sector coplanarity)
// ---------------------------------------------------------------------------

using hcg_tuple = std::tuple<std::vector<half_chain>, std::vector<std::size_t>,
                         std::vector<std::size_t>, std::vector<half_chain>,
                         std::vector<std::pair<std::size_t, std::size_t>>>;

hcg_tuple build_half_chain_graph(const chain_build_result& chains,
                                 const std::vector<point>& hot_pixels);

// ---------------------------------------------------------------------------
// Exterior face detection via ray casting
// ---------------------------------------------------------------------------

// Find exterior face info per connected component.
// Returns: the half-chain id in the exterior face (per component),
//   and all ray-coplanar pairs.
using fe_tuple = std::tuple<std::vector<std::size_t>,
                        std::vector<std::pair<std::size_t, std::size_t>>>;

fe_tuple find_exterior(const chain_build_result& chains,
                       const std::vector<point>& hot_pixels,
                       const std::vector<half_chain>& sorted_half_chains,
                       const std::vector<std::size_t>& half_chain_begin,
                       const std::vector<std::size_t>& half_chain_end);

// ---------------------------------------------------------------------------
// Winding numbers (coplanarity graph + DFS propagation from exterior seeds)
// Coplanarity sources: sector adjacency + ray.
// Propagation: for consecutive a,b at vertex CCW:
//   W(face(b)) = W(face(a)) + power(b)
// Returns: winding for each half-chain.
// ---------------------------------------------------------------------------

std::vector<int> compute_winding(
    const chain_build_result& chains,
    const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
    const std::vector<std::size_t>& exterior_half_chains);

// ---------------------------------------------------------------------------
// Build output polygons (filter → cancel → rebuild next → trace rings)
// ---------------------------------------------------------------------------
inline multi_polygon build_output(
    const chain_build_result& chains,
    const std::vector<point>& hot_pixels,
    std::vector<half_chain> next_half_chain,
    const std::vector<int>& winding,
    std::vector<std::pair<std::size_t, std::size_t>> exteral_coface_pairs,
    auto filter_fn)
{
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
    std::size_t num_chains = chains.offsets.size() - 1;

    // ----------------------------
    // 1. filter + dual cancellation
    // ----------------------------
    std::vector<bool> survive(num_half_chains);
    for (std::size_t i = 0; i < num_half_chains; i++)
        survive[i] = filter_fn(winding[i]);

    for (std::size_t i = 0; i < num_chains; i++) {
        std::size_t f = 2 * i, r = 2 * i + 1;
        if (survive[f] && survive[r]) {
            survive[f] = survive[r] = false;
            exteral_coface_pairs.emplace_back(f, r);
        }
    }

    // 修正 next_half_chain
    for (std::size_t i = 0; i < num_half_chains; i++) {
        if (!survive[i]) continue;
        while (!survive[next_half_chain[i].id])
            next_half_chain[i] =
                next_half_chain[next_half_chain[i].dual().id];
    }

    auto t1 = clock::now();

    // ----------------------------
    // 2. trace rings（核心变化）
    // ----------------------------
    std::vector<bool> visited(num_half_chains, false);
    std::vector<std::size_t> hc_to_ring(num_half_chains, (size_t)-1);
    std::vector<ring> rings;

    for (std::size_t i = 0; i < num_half_chains; i++) {

        if (!survive[i] || visited[i]) continue;

        std::size_t ring_id = rings.size();
        ring current_ring;

        std::size_t start = i;
        std::size_t cur = i;

        do {
            visited[cur] = true;
            hc_to_ring[cur] = ring_id;

            auto h = half_chain{cur};
            auto chain_idx = h.chain_id();
            auto begin = chains.offsets[chain_idx];
            auto end   = chains.offsets[chain_idx + 1];

            if (h.is_forward()) {
                for (std::size_t k = begin; k < end - 1; k++)
                    current_ring.push_back(hot_pixels[chains.indices[k]]);
            } else {
                for (std::size_t k = end - 1; k > begin; k--)
                    current_ring.push_back(hot_pixels[chains.indices[k]]);
            }

            cur = next_half_chain[cur].id;

        } while (cur != start);

        if (!current_ring.empty())
            current_ring.push_back(current_ring.front());

        rings.push_back(std::move(current_ring));
    }

    auto t2 = clock::now();

    // ----------------------------
    // 3. ring-level connected components
    // ----------------------------
    std::vector<std::pair<std::size_t, std::size_t>> ring_edges;
    ring_edges.reserve(exteral_coface_pairs.size());

    for (auto [a, b] : exteral_coface_pairs) {
        std::size_t ra = hc_to_ring[a];
        std::size_t rb = hc_to_ring[b];

        if (ra != (size_t)-1 && rb != (size_t)-1 && ra != rb)
            ring_edges.emplace_back(ra, rb);
    }

    auto component_id = connected_components(rings.size(), ring_edges);

    auto t3 = clock::now();

    // ----------------------------
    // 4. bucket: ring → polygon
    // ----------------------------
    std::vector<std::pair<std::size_t, std::size_t>> ring_to_poly;
    for (std::size_t i = 0; i < rings.size(); i++)
        ring_to_poly.emplace_back(component_id[i], i);

    std::size_t num_polygons = 0;
    for (auto c : component_id)
        num_polygons = std::max(num_polygons, c + 1);

    auto [poly_begin, poly_end, poly_data] = bucket_sort(
        ring_to_poly,
        num_polygons,
        [](const auto& p) { return p.first; },
        [](const auto& p) { return p.second; });

    auto t4 = clock::now();

    // ----------------------------
    // 5. build polygons
    // ----------------------------
    multi_polygon result;

    for (std::size_t pid = 0; pid < num_polygons; pid++) {

        auto begin = poly_begin[pid];
        auto end   = poly_end[pid];

        if (begin == end) continue;

        // 找 outer（最左点）
        std::size_t outer_idx = poly_data[begin];
        int32_t min_x = std::numeric_limits<int32_t>::max();

        for (std::size_t i = begin; i < end; i++) {
            auto rid = poly_data[i];

            for (auto const& p : rings[rid]) {
                int32_t x = bg::get<0>(p);
                if (x < min_x) {
                    min_x = x;
                    outer_idx = rid;
                }
            }
        }

        polygon poly;
        poly.outer() = std::move(rings[outer_idx]);

        for (std::size_t i = begin; i < end; i++) {
            auto rid = poly_data[i];
            if (rid == outer_idx) continue;
            if (!rings[rid].empty())
                poly.inners().push_back(std::move(rings[rid]));
        }

        result.push_back(std::move(poly));
    }

    auto t5 = clock::now();

    auto ms = [](auto d) {
        return std::chrono::duration<double, std::milli>(d).count();
    };

    std::fprintf(stderr,
        "  [output] filter=%.1fms trace=%.1fms cc=%.1fms bucket=%.1fms build=%.1fms (rings=%zu)\n",
        ms(t1 - t0),
        ms(t2 - t1),
        ms(t3 - t2),
        ms(t4 - t3),
        ms(t5 - t4),
        rings.size());

    return result;
}

// ---------------------------------------------------------------------------
// Main pipeline
// ---------------------------------------------------------------------------

inline auto run_pipeline(std::vector<point> points, std::vector<std::size_t> offsets, auto filter) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    auto [hot_pixels, sorted_edges] = construct_edges_with_power(points, offsets);
    auto t1 = clock::now();

    auto chains = build_chains(sorted_edges, hot_pixels.size());
    auto t2 = clock::now();

    auto [sorted_half_chains, half_chain_begin, half_chain_end, next_half_chain, coplanar] =
        build_half_chain_graph(chains, hot_pixels);
    auto t3 = clock::now();

    auto [exterior_half_chains, ray_pairs] =
        find_exterior(chains, hot_pixels, sorted_half_chains, half_chain_begin, half_chain_end);
    auto t4 = clock::now();

    auto winding = compute_winding(chains, coplanar, ray_pairs, exterior_half_chains);
    auto t5 = clock::now();

    auto result = build_output(chains, hot_pixels, std::move(next_half_chain),
                               winding, std::move(ray_pairs), filter);
    auto t6 = clock::now();

    auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };
    std::fprintf(stderr,
        "[pipeline] edges=%.1fms chains=%.1fms half_graph=%.1fms exterior=%.1fms "
        "winding=%.1fms output=%.1fms total=%.1fms\n",
        ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), ms(t4 - t3),
        ms(t5 - t4), ms(t6 - t5), ms(t6 - t0));
    return result;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

inline auto collect_segments(auto... mp_list) {
    std::vector<point> points;
    std::vector<std::size_t> offsets;

    auto add_ring = [&](const auto& ring) {
        if (ring.empty()) return;
        if (!points.empty())
            offsets.push_back(points.size() - 1);
        for (const auto& pt : ring)
            points.push_back(pt);
    };

    auto add_mp = [&](const auto& mp) {
        for (const auto& poly : mp) {
            add_ring(poly.outer());
            for (const auto& inner : poly.inners())
                add_ring(inner);
        }
    };

    (add_mp(mp_list), ...);

    return std::pair{std::move(points), std::move(offsets)};
}

inline auto add(const auto& polygons1, const auto& polygons2) {
    auto [points, offsets] = collect_segments(polygons1, polygons2);
    return run_pipeline(std::move(points), std::move(offsets), [](int w) { return w > 0; });
}

inline auto intersection(const auto& polygons1, const auto& polygons2) {
    auto [points, offsets] = collect_segments(polygons1, polygons2);
    return run_pipeline(std::move(points), std::move(offsets), [](int w) { return w > 1; });
}

inline auto self_or(auto r) {
    auto [points, offsets] = collect_segments(r);
    return run_pipeline(std::move(points), std::move(offsets), [](int w) { return w > 0; });
}

inline auto xor_(const auto& polygons1, const auto& polygons2) {
    auto [points, offsets] = collect_segments(polygons1, polygons2);
    return run_pipeline(std::move(points), std::move(offsets), [](int w) { return w == 1; });
}

inline auto difference(const auto& polygons1, const auto& polygons2) {
    auto [points, offsets] = collect_segments(polygons1);

    auto add_mp_rev = [&](const auto& mp) {
        for (const auto& poly : mp) {
            auto add_ring_rev = [&](const auto& ring) {
                if (ring.empty()) return;
                if (!points.empty())
                    offsets.push_back(points.size() - 1);
                points.push_back(ring[0]);
                for (int i = (int)ring.size() - 2; i >= 0; i--)
                    points.push_back(ring[i]);
            };
            add_ring_rev(poly.outer());
            for (const auto& inner : poly.inners())
                add_ring_rev(inner);
        }
    };
    add_mp_rev(polygons2);

    return run_pipeline(std::move(points), std::move(offsets), [](int w) { return w > 0; });
}

} // namespace best_clipper
