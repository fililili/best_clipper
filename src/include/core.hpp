#pragma once

#include <algorithm>
#include <boost/container/flat_map.hpp>
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

inline multi_polygon build_output(const chain_build_result& chains,
                                   const std::vector<point>& hot_pixels,
                                   std::vector<half_chain> next_half_chain,
                                   const std::vector<int>& winding,
                                   const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
                                   const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
                                   auto filter_fn) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();
    std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;
    std::size_t num_chains = chains.offsets.size() - 1;

    // Filter half-chains by winding number
    std::vector<bool> survive(num_half_chains);
    for (std::size_t i = 0; i < num_half_chains; i++)
        survive[i] = filter_fn(winding[i]);

    // Dual cancellation: both fwd and rev survive → both dead
    std::vector<bool> dead(num_half_chains, false);
    for (std::size_t chain_idx = 0; chain_idx < num_chains; chain_idx++) {
        std::size_t forward_id = 2 * chain_idx, reverse_id = 2 * chain_idx + 1;
        if (survive[forward_id] && survive[reverse_id])
            dead[forward_id] = dead[reverse_id] = true;
    }

    // Update next_half_chain via indirection: for each surviving half-chain whose next
    // points to a dead half-chain, follow next_half_chain[dead.dual()] until a live
    // half-chain is found. Terminates because at least one half-chain per vertex cycle
    // survives (not all faces can be cancelled).
    for (std::size_t i = 0; i < num_half_chains; i++) {
        if (!survive[i] || dead[i]) continue;
        while (dead[next_half_chain[i].id])
            next_half_chain[i] = next_half_chain[next_half_chain[i].dual().id];
    }
    auto t1 = clock::now();

    // Build connected components: coplanar + ray + dual cancellation → same face
    std::vector<std::pair<std::size_t, std::size_t>> face_edges;
    face_edges.reserve(coplanar_pairs.size() + ray_pairs.size() + num_chains);
    face_edges.insert(face_edges.end(), coplanar_pairs.begin(), coplanar_pairs.end());
    face_edges.insert(face_edges.end(), ray_pairs.begin(), ray_pairs.end());
    for (std::size_t chain_idx = 0; chain_idx < num_chains; chain_idx++) {
        std::size_t forward_id = 2 * chain_idx, reverse_id = 2 * chain_idx + 1;
        if (dead[forward_id]) face_edges.emplace_back(forward_id, reverse_id);
    }
    auto component_id = connected_components(num_half_chains, face_edges);
    auto t2 = clock::now();

    // Group surviving non-dead half-chains by face
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
    auto t3 = clock::now();

    // For each face, trace its half-chains into a polygon
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
            int32_t min_x = std::numeric_limits<int32_t>::max();
            for (std::size_t ring_index = 0; ring_index < rings.size(); ring_index++) {
                for (const auto& point_value : rings[ring_index]) {
                    int32_t x = bg::get<0>(point_value);
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

    auto t4 = clock::now();
    auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };
    std::fprintf(stderr, "  [output] filter=%.1fms cc=%.1fms bucket=%.1fms trace=%.1fms (faces=%zu)\n",
        ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), ms(t4 - t3), num_faces);
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
                               winding, coplanar, ray_pairs, filter);
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
