#pragma once

// Internal header: implementation details only used by core.cpp.
// Must be included AFTER core.hpp (depends on point, segment, int128_t, etc.).

#include <cstdint>
#include <optional>
#include <utility>

namespace best_clipper {

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

struct edge_t { std::size_t start; std::size_t end; };

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

inline bool product_gt(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
    return a * b > c * d;
}

inline bool less_by_direction(point source, point target1, point target2) {
    int64_t dx1 = (int64_t)bg::get<0>(target1) - (int64_t)bg::get<0>(source);
    int64_t dy1 = (int64_t)bg::get<1>(target1) - (int64_t)bg::get<1>(source);
    int64_t dx2 = (int64_t)bg::get<0>(target2) - (int64_t)bg::get<0>(source);
    int64_t dy2 = (int64_t)bg::get<1>(target2) - (int64_t)bg::get<1>(source);

    auto quadrant = [](int64_t dx, int64_t dy) {
        if (dx > 0 && dy >= 0) return 0;
        if (dx <= 0 && dy > 0) return 1;
        if (dx < 0 && dy <= 0) return 2;
        return 3;
    };
    int q1 = quadrant(dx1, dy1), q2 = quadrant(dx2, dy2);
    if (q1 != q2) return q1 < q2;

    uint64_t ux1 = dx1 >= 0 ? (uint64_t)dx1 : (uint64_t)(-dx1);
    uint64_t uy1 = dy1 >= 0 ? (uint64_t)dy1 : (uint64_t)(-dy1);
    uint64_t ux2 = dx2 >= 0 ? (uint64_t)dx2 : (uint64_t)(-dx2);
    uint64_t uy2 = dy2 >= 0 ? (uint64_t)dy2 : (uint64_t)(-dy2);
    bool cross_positive = (q1 % 2 == 0) ? product_gt(ux1, uy2, uy1, ux2)
                                        : product_gt(uy1, ux2, ux1, uy2);
    return cross_positive;
}

struct less_by_segment {
    int32_t dx, dy;

    less_by_segment(const segment& s)
        : dx(bg::get<1, 0>(s) - bg::get<0, 0>(s))
        , dy(bg::get<1, 1>(s) - bg::get<0, 1>(s)) {}

    bool operator()(point p1, point p2) const {
        auto x1 = bg::get<0>(p1), x2 = bg::get<0>(p2);
        if (dx != 0) {
            if (x1 != x2) return dx > 0 ? x1 < x2 : x1 > x2;
            auto y1 = bg::get<1>(p1), y2 = bg::get<1>(p2);
            return dy >= 0 ? y1 < y2 : y1 > y2;
        }
        auto y1 = bg::get<1>(p1), y2 = bg::get<1>(p2);
        return dy > 0 ? y1 < y2 : y1 > y2;
    }
};

inline int128_t cross_i128(int64_t dx1, int64_t dy1, int64_t dx2, int64_t dy2) {
    return (int128_t)dx1 * dy2 - (int128_t)dy1 * dx2;
}

inline bool cross_eq(int64_t dx1, int64_t dy1, int64_t dx2, int64_t dy2) {
    if ((dx1 >= 0) != (dy2 >= 0)) return false;
    if ((dy1 >= 0) != (dx2 >= 0)) return false;
    uint64_t a = (uint64_t)(dx1 >= 0 ? dx1 : -dx1) * (uint64_t)(dy2 >= 0 ? dy2 : -dy2);
    uint64_t b = (uint64_t)(dy1 >= 0 ? dy1 : -dy1) * (uint64_t)(dx2 >= 0 ? dx2 : -dx2);
    return a == b;
}

inline std::optional<point> get_intersection(segment s1, segment s2) {
    int64_t x1 = bg::get<0, 0>(s1), y1 = bg::get<0, 1>(s1);
    int64_t x2 = bg::get<1, 0>(s1), y2 = bg::get<1, 1>(s1);
    int64_t x3 = bg::get<0, 0>(s2), y3 = bg::get<0, 1>(s2);
    int64_t x4 = bg::get<1, 0>(s2), y4 = bg::get<1, 1>(s2);
    int64_t dx1 = x2 - x1, dy1 = y2 - y1;
    int64_t dx2 = x4 - x3, dy2 = y4 - y3;

    // Fast path: axis-aligned segments — no int128 needed.
    if (dx1 == 0 && dy2 == 0) {
        // s1 vertical, s2 horizontal → intersection at (x1, y3)
        if (x1 <= std::min(x3, x4) || x1 >= std::max(x3, x4)) return {};
        if (y3 <= std::min(y1, y2) || y3 >= std::max(y1, y2)) return {};
        return point{(int32_t)x1, (int32_t)y3};
    }
    if (dy1 == 0 && dx2 == 0) {
        // s1 horizontal, s2 vertical → intersection at (x3, y1)
        if (x3 <= std::min(x1, x2) || x3 >= std::max(x1, x2)) return {};
        if (y1 <= std::min(y3, y4) || y1 >= std::max(y3, y4)) return {};
        return point{(int32_t)x3, (int32_t)y1};
    }

    if (cross_eq(dx1, dy1, dx2, dy2)) return {};

    int128_t d = cross_i128(dx1, dy1, dx2, dy2);
    int128_t t1_num = cross_i128(x3 - x1, y3 - y1, dx2, dy2);
    int128_t t2_num = cross_i128(x3 - x1, y3 - y1, dx1, dy1);
    if (d < 0) { d = -d; t1_num = -t1_num; t2_num = -t2_num; }
    if (t1_num <= 0 || t1_num >= d || t2_num <= 0 || t2_num >= d) return {};

    auto snap = [&](int64_t x, int64_t dx, int128_t num, int128_t den) -> int32_t {
        int128_t num_dx = num * dx;
        int128_t n = 2 * num_dx + den - 1;
        int128_t dd = 2 * den;
        int64_t q = (int64_t)(n >= 0 ? n / dd : (n - dd + 1) / dd);
        return (int32_t)(x + q);
    };
    return point{snap(x1, dx1, t1_num, d), snap(y1, dy1, t1_num, d)};
}

inline bool is_point_on_segment(point p, segment s) {
    int64_t x = bg::get<0>(p), y = bg::get<1>(p);
    int64_t x1 = bg::get<0, 0>(s), y1 = bg::get<0, 1>(s);
    int64_t x2 = bg::get<1, 0>(s), y2 = bg::get<1, 1>(s);
    int64_t dxs = x2 - x1, dys = y2 - y1;
    {
        int64_t dxp = x - x1, dyp = y - y1;
        int sign_a = (dxp ^ dxs) >= 0 ? 1 : -1;
        int sign_b = (dyp ^ dys) >= 0 ? 1 : -1;
        uint64_t abs_a = (uint64_t)(dxp >= 0 ? dxp : -dxp) * (uint64_t)(dxs >= 0 ? dxs : -dxs);
        uint64_t abs_b = (uint64_t)(dyp >= 0 ? dyp : -dyp) * (uint64_t)(dys >= 0 ? dys : -dys);
        bool dot_lt_0;
        if (sign_a == sign_b) dot_lt_0 = sign_a < 0 && (abs_a > 0 || abs_b > 0);
        else dot_lt_0 = sign_a > 0 ? abs_b > abs_a : abs_a > abs_b;
        if (dot_lt_0) return false;
    }
    {
        int64_t dxp = x2 - x, dyp = y2 - y;
        int sign_a = (dxp ^ dxs) >= 0 ? 1 : -1;
        int sign_b = (dyp ^ dys) >= 0 ? 1 : -1;
        uint64_t abs_a = (uint64_t)(dxp >= 0 ? dxp : -dxp) * (uint64_t)(dxs >= 0 ? dxs : -dxs);
        uint64_t abs_b = (uint64_t)(dyp >= 0 ? dyp : -dyp) * (uint64_t)(dys >= 0 ? dys : -dys);
        bool dot_lt_0;
        if (sign_a == sign_b) dot_lt_0 = sign_a < 0 && (abs_a > 0 || abs_b > 0);
        else dot_lt_0 = sign_a > 0 ? abs_b > abs_a : abs_a > abs_b;
        if (dot_lt_0) return false;
    }
    int64_t A = 2 * x + 1 - x1 - x2, B = 2 * x - 1 - x1 - x2;
    int64_t C = 2 * y + 1 - y1 - y2, D = 2 * y - 1 - y1 - y2;
    auto cmp_gt = [](int64_t a, int64_t b, int64_t c, int64_t d) {
        int sign_a = (a ^ b) >= 0 ? 1 : -1;
        int sign_b = (c ^ d) >= 0 ? 1 : -1;
        if (sign_a != sign_b) return sign_a > sign_b;
        uint64_t abs_a = (uint64_t)(a >= 0 ? a : -a) * (uint64_t)(b >= 0 ? b : -b);
        uint64_t abs_b = (uint64_t)(c >= 0 ? c : -c) * (uint64_t)(d >= 0 ? d : -d);
        return sign_a > 0 ? abs_a > abs_b : abs_a < abs_b;
    };
    auto cmp_ge = [&](int64_t a, int64_t b, int64_t c, int64_t d) {
        return !cmp_gt(c, d, a, b);
    };
    auto cmp_lt = [&](int64_t a, int64_t b, int64_t c, int64_t d) {
        return cmp_gt(c, d, a, b);
    };
    auto cmp_le = [&](int64_t a, int64_t b, int64_t c, int64_t d) {
        return !cmp_gt(a, b, c, d);
    };
    if (cmp_gt(A, dys, C, dxs) && cmp_ge(A, dys, D, dxs) && cmp_ge(B, dys, C, dxs) && cmp_ge(B, dys, D, dxs)) return false;
    if (cmp_lt(A, dys, C, dxs) && cmp_le(A, dys, D, dxs) && cmp_le(B, dys, C, dxs) && cmp_le(B, dys, D, dxs)) return false;
    return true;
}

} // namespace best_clipper
