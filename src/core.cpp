#include "include/core.hpp"
#include "core_detail.hpp"
#include "include/uniform_grid.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>
#include <vector>

namespace best_clipper {

namespace bg = boost::geometry;

// ---------------------------------------------------------------------------
// connected_components
// ---------------------------------------------------------------------------

std::vector<std::size_t> connected_components(
    std::size_t n,
    const std::vector<std::pair<std::size_t, std::size_t>>& edges) {

    std::vector<std::size_t> parent(n);
    for (std::size_t i = 0; i < n; i++) parent[i] = i;

    auto find = [&](std::size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    for (auto [a, b] : edges) {
        auto ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    }

    std::vector<std::size_t> comp(n);
    std::vector<std::size_t> root_to_id(n, ~0ULL);
    std::size_t num_components = 0;
    for (std::size_t v = 0; v < n; v++) {
        auto r = find(v);
        if (root_to_id[r] == ~0ULL) root_to_id[r] = num_components++;
        comp[v] = root_to_id[r];
    }

    return comp;
}

// ---------------------------------------------------------------------------
// construct_graph
// ---------------------------------------------------------------------------

std::tuple<std::vector<edge_t>, std::vector<point>>
construct_graph(const std::vector<point>& points, const std::vector<std::size_t>& offsets) {
    // Build grid: point index i is the segment (points[i], points[i+1]).
    // Offset segments get an invalid bbox so the grid skips them naturally.
    std::vector<box> seg_boxes;
    std::vector<std::size_t> seg_idx;
    seg_boxes.reserve(points.size() - 1);
    seg_idx.reserve(points.size() - 1);
    {
        std::size_t start = 0;
        for (auto off : offsets) {
            for (std::size_t i = start; i < off; i++) {
                seg_boxes.push_back(bg::return_envelope<box>(segment{points[i], points[i + 1]}));
                seg_idx.push_back(i);
            }
            seg_boxes.push_back(box{point{1, 1}, point{0, 0}});
            seg_idx.push_back(off);
            start = off + 1;
        }
        for (std::size_t i = start; i + 1 < points.size(); i++) {
            seg_boxes.push_back(bg::return_envelope<box>(segment{points[i], points[i + 1]}));
            seg_idx.push_back(i);
        }
    }
    best_clipper::uniform_grid::grid segments_box_grid(seg_boxes);

    // Hot pixels from all points
    std::vector<point> hot_pixels = points;

    // Find intersections between segments (skip offset boundaries)
    {
        std::size_t start = 0;
        for (auto off : offsets) {
            for (std::size_t i = start; i < off; i++) {
                auto si = segment{points[i], points[i + 1]};
                box box_i = bg::return_envelope<box>(si);
                segments_box_grid.query_intersects(box_i, [&](std::size_t pos) {
                    std::size_t j = seg_idx[pos];
                    if (i <= j + 1) return; // adjacent edges share a vertex, cannot produce a new intersection
                    if (auto p = get_intersection(si, segment{points[j], points[j + 1]}))
                        hot_pixels.push_back(p.value());
                });
            }
            start = off + 1;
        }
        for (std::size_t i = start; i + 1 < points.size(); i++) {
            auto si = segment{points[i], points[i + 1]};
            box box_i = bg::return_envelope<box>(si);
            segments_box_grid.query_intersects(box_i, [&](std::size_t pos) {
                std::size_t j = seg_idx[pos];
                if (i <= j + 1) return;
                if (auto p = get_intersection(si, segment{points[j], points[j + 1]}))
                    hot_pixels.push_back(p.value());
            });
        }
    }

    // Dedup hot_pixels: group by grid cell, sort within cell, then unique.
    // Better cache locality than global sort since spatially-close points
    // (which are more likely to be duplicates) fall in the same cell.
    {
        auto& g = segments_box_grid;
        std::size_t num_cells = (std::size_t)(g._x_cells * g._y_cells);
        std::vector<std::size_t> pixel_cell(hot_pixels.size());
        for (std::size_t pi = 0; pi < hot_pixels.size(); pi++) {
            auto& p = hot_pixels[pi];
            pixel_cell[pi] = (std::size_t)(((bg::get<1>(p) - g._min_y) / g._cell_size) * g._x_cells
                                          + (bg::get<0>(p) - g._min_x) / g._cell_size);
        }
        std::vector<std::size_t> cell_counts(num_cells, 0);
        for (auto c : pixel_cell) cell_counts[c]++;
        std::vector<std::size_t> begins(num_cells + 1, 0);
        for (std::size_t c = 0; c < num_cells; c++)
            begins[c + 1] = begins[c] + cell_counts[c];
        std::vector<std::size_t> cell_items(hot_pixels.size());
        auto cursors = begins;
        for (std::size_t pi = 0; pi < hot_pixels.size(); pi++)
            cell_items[cursors[pixel_cell[pi]]++] = pi;
        std::vector<point> unique;
        unique.reserve(hot_pixels.size());
        for (std::size_t c = 0; c < num_cells; c++) {
            auto b = begins[c], e = begins[c + 1];
            std::sort(cell_items.begin() + b, cell_items.begin() + e,
                [&](std::size_t a, std::size_t b) {
                    return std::pair{bg::get<0>(hot_pixels[a]), bg::get<1>(hot_pixels[a])} <
                           std::pair{bg::get<0>(hot_pixels[b]), bg::get<1>(hot_pixels[b])};
                });
            for (auto it = cell_items.begin() + b; it != cell_items.begin() + e; ++it) {
                if (it == cell_items.begin() + b ||
                    !bg::equals(hot_pixels[*it], hot_pixels[*(it - 1)]))
                    unique.push_back(hot_pixels[*it]);
            }
        }
        hot_pixels = std::move(unique);
    }

    // Build segment-pixel pairs
    std::vector<std::pair<std::size_t, std::size_t>> segment_pixel_pairs;
    std::vector<std::size_t> last_seen(points.size(), ~0ULL);
    for (std::size_t pi = 0; pi < hot_pixels.size(); pi++) {
        int32_t x = bg::get<0>(hot_pixels[pi]), y = bg::get<1>(hot_pixels[pi]);
        auto min_corner = point{x > INT32_MIN ? x - 1 : INT32_MIN, y > INT32_MIN ? y - 1 : INT32_MIN};
        auto max_corner = point{x < INT32_MAX ? x + 1 : INT32_MAX, y < INT32_MAX ? y + 1 : INT32_MAX};
        segments_box_grid.query_intersects(box{min_corner, max_corner}, [&](std::size_t pos) {
            std::size_t seg_start = seg_idx[pos];
            if (last_seen[seg_start] == pi) return;
            last_seen[seg_start] = pi;
            if (is_point_on_segment(hot_pixels[pi], segment{points[seg_start], points[seg_start + 1]}))
                segment_pixel_pairs.emplace_back(seg_start, pi);
        });
    }

    // Build edges: sort pixels along each segment, connect adjacent pairs
    std::vector<edge_t> edges;
    {
        auto [segments_begin, segments_end, pixels] = bucket_sort(
            segment_pixel_pairs, points.size(), [](auto val) { return val.first; },
            [](auto val) { return val.second; });
        std::size_t start = 0;
        for (auto off : offsets) {
            for (std::size_t i = start; i < off; i++) {
                auto pixel_begin = std::begin(pixels) + segments_begin[i],
                     pixel_end   = std::begin(pixels) + segments_end[i];
                auto seg = segment{points[i], points[i + 1]};
                std::sort(pixel_begin, pixel_end, [&](auto pi, auto pj) {
                    return less_by_segment(seg)(hot_pixels[pi], hot_pixels[pj]);
                });
                for (; pixel_begin != pixel_end - 1; pixel_begin++)
                    edges.emplace_back(*pixel_begin, *std::next(pixel_begin));
            }
            start = off + 1;
        }
        for (std::size_t i = start; i + 1 < points.size(); i++) {
            auto pixel_begin = std::begin(pixels) + segments_begin[i],
                 pixel_end   = std::begin(pixels) + segments_end[i];
            auto seg = segment{points[i], points[i + 1]};
            std::sort(pixel_begin, pixel_end, [&](auto pi, auto pj) {
                return less_by_segment(seg)(hot_pixels[pi], hot_pixels[pj]);
            });
            for (; pixel_begin != pixel_end - 1; pixel_begin++)
                edges.emplace_back(*pixel_begin, *std::next(pixel_begin));
        }
    }
    return {std::move(edges), std::move(hot_pixels)};
}

// ---------------------------------------------------------------------------
// edges_to_power
// ---------------------------------------------------------------------------

std::vector<edge_with_power_t> edges_to_power(std::vector<edge_t> edges) {
    std::vector<edge_with_power_t> r(edges.size());
    for (std::size_t i = 0; i < r.size(); i++) {
        auto s = edges[i].start, e = edges[i].end;
        if (s < e) r[i] = {s, e, 1};
        else      r[i] = {e, s, -1};
    }
    return r;
}

// ---------------------------------------------------------------------------
// unique_edges
// ---------------------------------------------------------------------------

std::vector<edge_with_power_t> unique_edges(std::vector<edge_with_power_t> edges,
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

// ---------------------------------------------------------------------------
// construct_edges_with_power
// ---------------------------------------------------------------------------

std::tuple<std::vector<point>, std::vector<edge_with_power_t>>
construct_edges_with_power(const std::vector<point>& points, const std::vector<std::size_t>& offsets) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto [edges, hot_pixels] = construct_graph(points, offsets);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto r = std::tuple{std::move(hot_pixels),
                        unique_edges(edges_to_power(std::move(edges)), hot_pixels.size())};
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };
    std::fprintf(stderr, "  [edges] graph=%.1fms dedup=%.1fms (hp=%zu)\n",
        ms(t1 - t0), ms(t2 - t1), std::get<0>(r).size());
    return r;
}

// ---------------------------------------------------------------------------
// build_chains
// ---------------------------------------------------------------------------

chain_build_result build_chains(const std::vector<edge_with_power_t>& sorted_edges,
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
// build_half_chain_graph
// ---------------------------------------------------------------------------

using hcg_tuple = std::tuple<std::vector<half_chain>, std::vector<std::size_t>,
                             std::vector<std::size_t>, std::vector<half_chain>,
                             std::vector<std::pair<std::size_t, std::size_t>>>;

hcg_tuple build_half_chain_graph(const chain_build_result& chains,
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

    return hcg_tuple{
        std::move(sorted_half_chains),
        std::vector<std::size_t>(begin_loc.begin(), begin_loc.end()),
        std::vector<std::size_t>(end_loc.begin(), end_loc.end()),
        std::move(next_half_chain), std::move(coplanar)};
}

// ---------------------------------------------------------------------------
// cast_ray_minus_x
// ---------------------------------------------------------------------------

struct chain_seg {
    int64_t x1, y1, dx, dy;
    std::size_t fwd_hc, rev_hc;
};

std::vector<std::pair<std::size_t, std::size_t>> cast_ray_minus_x(
    std::size_t v,
    const std::vector<point>& hot_pixels,
    const chain_build_result& chains,
    const std::vector<half_chain>& sorted_half_chains,
    const std::vector<std::size_t>& half_chain_begin,
    const std::vector<std::size_t>& half_chain_end,
    const best_clipper::uniform_grid::grid& seg_grid,
    const std::vector<chain_seg>& seg_data) {

    auto range_begin = half_chain_begin[v], range_end = half_chain_end[v];
    if (range_begin == range_end) return {};

    // Find the half-chain at v whose right face contains the -x ray direction.
    point vertex_point = hot_pixels[v];
    int32_t vx = bg::get<0>(vertex_point);
    int32_t vy = bg::get<1>(vertex_point);
    point ray_point{vx > INT32_MIN ? vx - 1 : INT32_MIN, vy};
    half_chain vertex_half_chain;
    bool found = false;
    for (auto it = range_begin; it < range_end; it++) {
        auto h = sorted_half_chains[it];
        auto half_chain_point = hot_pixels[h.next_along_source(chains)];
        if (less_by_direction(vertex_point, ray_point, half_chain_point)) {
            vertex_half_chain = h;
            found = true;
            break;
        }
    }
    if (!found) vertex_half_chain = sorted_half_chains[range_begin];

    // Cast ray in -x direction, find the nearest intersected half-chain.
    int64_t ray_y = vy;
    int64_t min_x = vx;
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
            int sign_a = (dx ^ hit_dy) >= 0 ? 1 : -1;
            int sign_b = (best_dx ^ dy) >= 0 ? 1 : -1;
            bool same_sign = sign_a == sign_b;
            uint64_t abs_a = (uint64_t)(dx >= 0 ? dx : -dx) * (uint64_t)(hit_dy >= 0 ? hit_dy : -hit_dy);
            uint64_t abs_b = (uint64_t)(best_dx >= 0 ? best_dx : -best_dx) * (uint64_t)(dy >= 0 ? dy : -dy);
            bool slope_gt;
            if (!same_sign) slope_gt = sign_a > sign_b;
            else if (sign_a > 0) slope_gt = abs_a > abs_b;
            else slope_gt = abs_a < abs_b;
            if (slope_gt) {
                better = true;
            } else if (same_sign && abs_a == abs_b) {
                uint64_t da = ray_y - y_low;
                uint64_t db = ray_y - hit_y_low;
                uint64_t abs_dy = dy > 0 ? (uint64_t)dy : (uint64_t)(-dy);
                uint64_t abs_hdy = hit_dy > 0 ? (uint64_t)hit_dy : (uint64_t)(-hit_dy);
                if (da * abs_hdy < db * abs_dy) {
                    better = true;
                } else if (da * abs_hdy == db * abs_dy && abs_dy > abs_hdy) {
                    better = true;
                }
            }
        }
        if (better) {
            best_x = ix; best_dx = dx; hit_dy = dy; hit_id = half_chain_id;
            hit_y_low = y_low;
        }
    };

    auto intersect_x = [&](int64_t x1, int64_t y1, int64_t dx, int64_t dy) -> int64_t {
        return (int64_t)(((int128_t)x1 * dy + (int128_t)(ray_y - y1) * dx) / dy);
    };

    auto query_box = box{point{INT32_MIN, (int32_t)ray_y}, point{vx, (int32_t)ray_y}};
    seg_grid.query_intersects(query_box, [&](std::size_t idx) {
        auto& seg = seg_data[idx];
        int64_t y2 = seg.y1 + seg.dy;
        if (seg.y1 <= ray_y && ray_y < y2) {
            int64_t ix = intersect_x(seg.x1, seg.y1, seg.dx, seg.dy);
            try_edge(ix, seg.dx, seg.dy, seg.y1, seg.fwd_hc);
        } else if (y2 <= ray_y && ray_y < seg.y1) {
            int64_t ix = intersect_x(seg.x1, seg.y1, seg.dx, seg.dy);
            try_edge(ix, -seg.dx, -seg.dy, y2, seg.rev_hc);
        }
    });

    // Build the coplanar pair (RIGHT-face convention).
    std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;
    if (hit_id != ~0ULL) {
        std::size_t hit_side = (hit_dy > 0) ? hit_id : (hit_id ^ 1);
        ray_pairs.emplace_back(vertex_half_chain.id, hit_side);
    } else {
        ray_pairs.emplace_back(vertex_half_chain.id, vertex_half_chain.id);
    }
    return ray_pairs;
}

// ---------------------------------------------------------------------------
// find_exterior
// ---------------------------------------------------------------------------

using fe_tuple = std::tuple<std::vector<std::size_t>,
                            std::vector<std::pair<std::size_t, std::size_t>>>;

fe_tuple find_exterior(const chain_build_result& chains,
                       const std::vector<point>& hot_pixels,
                       const std::vector<half_chain>& sorted_half_chains,
                       const std::vector<std::size_t>& half_chain_begin,
                       const std::vector<std::size_t>& half_chain_end) {
    std::size_t num_vertices = hot_pixels.size();

    std::vector<std::pair<std::size_t, std::size_t>> vertex_edges;
    for (auto h : sorted_half_chains)
        vertex_edges.emplace_back(h.source_node(chains), h.target_node(chains));
    for (std::size_t c = 0; c + 1 < chains.offsets.size(); c++) {
        auto chain_begin_idx = chains.offsets[c], chain_end_idx = chains.offsets[c + 1];
        for (std::size_t k = chain_begin_idx; k + 1 < chain_end_idx; k++)
            vertex_edges.emplace_back(chains.indices[k], chains.indices[k + 1]);
    }

    auto t_cc0 = std::chrono::high_resolution_clock::now();
    auto component_id = connected_components(num_vertices, vertex_edges);
    auto t_cc1 = std::chrono::high_resolution_clock::now();

    std::size_t num_components = 0;
    for (auto c : component_id) num_components = std::max(num_components, c + 1);
    std::vector<std::vector<std::size_t>> vertex_components(num_components);
    for (std::size_t v = 0; v < num_vertices; v++)
        vertex_components[component_id[v]].push_back(v);
    auto t_group = std::chrono::high_resolution_clock::now();
    auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };
    std::fprintf(stderr, "  [exterior] cc=%.1fms group=%.1fms (n=%zu m=%zu comps=%zu)\n",
        ms(t_cc1 - t_cc0), ms(t_group - t_cc1), num_vertices, vertex_edges.size(), num_components);

    auto t_data0 = std::chrono::high_resolution_clock::now();
    std::vector<chain_seg> seg_data;
    std::vector<box> seg_boxes;
    size_t total_edges = 0;
    for (size_t ci = 0; ci + 1 < chains.offsets.size(); ci++)
        total_edges += chains.offsets[ci + 1] - chains.offsets[ci] - 1;
    seg_data.reserve(total_edges);
    seg_boxes.reserve(total_edges);
    for (size_t ci = 0; ci + 1 < chains.offsets.size(); ci++) {
        auto cbi = chains.offsets[ci], cei = chains.offsets[ci + 1];
        for (size_t k = cbi; k + 1 < cei; k++) {
            point p1 = hot_pixels[chains.indices[k]];
            point p2 = hot_pixels[chains.indices[k + 1]];
            int32_t x1 = bg::get<0>(p1), y1 = bg::get<1>(p1);
            int32_t x2 = bg::get<0>(p2), y2 = bg::get<1>(p2);
            seg_data.push_back({(int64_t)x1, (int64_t)y1,
                                (int64_t)x2 - (int64_t)x1,
                                (int64_t)y2 - (int64_t)y1,
                                2 * ci, 2 * ci + 1});
            seg_boxes.push_back(box{point{std::min(x1, x2), std::min(y1, y2)},
                                    point{std::max(x1, x2), std::max(y1, y2)}});
        }
    }
    best_clipper::uniform_grid::grid seg_grid(seg_boxes);
    auto t_data1 = std::chrono::high_resolution_clock::now();

    std::vector<std::size_t> exterior_half_chains;
    std::vector<std::pair<std::size_t, std::size_t>> ray_pairs;

    for (auto& vertex_component : vertex_components) {
        std::size_t leftmost_vertex = ~0ULL;
        int32_t min_x = std::numeric_limits<int32_t>::max();
        for (auto vertex : vertex_component) {
            if (half_chain_begin[vertex] == half_chain_end[vertex]) continue;
            int32_t x = bg::get<0>(hot_pixels[vertex]);
            if (x < min_x) { min_x = x; leftmost_vertex = vertex; }
        }
        if (leftmost_vertex == ~0ULL) continue;

        auto ray_pairs_result = cast_ray_minus_x(leftmost_vertex, hot_pixels, chains, sorted_half_chains, half_chain_begin, half_chain_end, seg_grid, seg_data);
        for (auto& p : ray_pairs_result) {
            if (p.first == p.second) {
                exterior_half_chains.push_back(p.first);
            }
            ray_pairs.push_back(std::move(p));
        }
    }

    auto t_rays = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "  [exterior] data=%.1fms grid=%.1fms rays=%.1fms\n",
        ms(t_data0 - t_group), ms(t_data1 - t_data0), ms(t_rays - t_data1));
    return fe_tuple{std::move(exterior_half_chains), std::move(ray_pairs)};
}

// ---------------------------------------------------------------------------
// compute_winding
// ---------------------------------------------------------------------------

std::vector<int> compute_winding(
    const chain_build_result& chains,
    const std::vector<std::pair<std::size_t, std::size_t>>& coplanar_pairs,
    const std::vector<std::pair<std::size_t, std::size_t>>& ray_pairs,
    const std::vector<std::size_t>& exterior_half_chains) {

    auto t0 = std::chrono::high_resolution_clock::now();
    std::size_t num_half_chains = (chains.offsets.size() - 1) * 2;

    // Build adjacency for coplanar + ray edges only.
    // Dual edges (i ↔ i^1) handled explicitly in DFS to avoid
    // bloating the adjacency by 2×num_half_chains entries.
    std::vector<edge_with_power_t> edges;
    edges.reserve((coplanar_pairs.size() + ray_pairs.size()) * 2);
    for (auto [a, b] : coplanar_pairs) {
        edges.emplace_back(a, b, 0);
        edges.emplace_back(b, a, 0);
    }
    for (auto [a, b] : ray_pairs) {
        edges.emplace_back(a, b, 0);
        edges.emplace_back(b, a, 0);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    auto [adj_begin, adj_end, adjacency] = bucket_sort(
        edges, num_half_chains,
        [](const edge_with_power_t& e) { return e.start; },
        [](const edge_with_power_t& e) { return std::pair{e.end, e.power}; });
    auto t2 = std::chrono::high_resolution_clock::now();

    constexpr int UNKNOWN = std::numeric_limits<int>::max() / 2;
    std::vector<int> winding(num_half_chains, UNKNOWN);
    for (auto exterior : exterior_half_chains) winding[exterior] = 0;

    std::vector<std::size_t> stack(exterior_half_chains.begin(), exterior_half_chains.end());
    while (!stack.empty()) {
        auto u = stack.back(); stack.pop_back();
        for (auto j = adj_begin[u]; j < adj_end[u]; j++) {
            auto [v, diff] = adjacency[j];
            if (winding[v] == UNKNOWN) {
                winding[v] = winding[u] + diff;
                stack.push_back(v);
            }
        }
        // Dual edge: i ↔ i^1 with winding difference = -power of i
        auto dual = u ^ 1;
        if (winding[dual] == UNKNOWN) {
            winding[dual] = winding[u] - half_chain{u}.power(chains);
            stack.push_back(dual);
        }
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };
    std::fprintf(stderr, "  [winding] build=%.1fms sort=%.1fms dfs=%.1fms (hc=%zu)\n",
        ms(t1 - t0), ms(t2 - t1), ms(t3 - t2), num_half_chains);
    return winding;
}

} // namespace best_clipper
