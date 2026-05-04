#pragma once
/// Boolean operations with automatic int32 ↔ uint32 conversion.
///
/// int32_t coordinates are biased by XOR with 0x80000000 (equivalent to +2^31)
/// before entering the uint32_t algorithm, preserving total ordering across the
/// signed/unsigned boundary. Results are unbiased on return.

#include <boost/geometry.hpp>
#include <cstdint>

#include "core.hpp"

namespace bg = boost::geometry;

using point_s32 = bg::model::d2::point_xy<int32_t>;
using segment_s32 = bg::model::segment<point_s32>;
using box_s32 = bg::model::box<point_s32>;
using ring_s32 = bg::model::ring<point_s32>;
using polygon_s32 = bg::model::polygon<point_s32>;
using multi_polygon_s32 = bg::model::multi_polygon<polygon_s32>;
using multi_polygon_u32 = best_clipper::multi_polygon;

using best_clipper::point;
using best_clipper::ring;
using best_clipper::polygon;

constexpr uint32_t BIAS = 0x80000000;

inline uint32_t bias_coord(int32_t v) { return (uint32_t)v ^ BIAS; }
inline int32_t unbias_coord(uint32_t v) { return (int32_t)(v ^ BIAS); }

inline point bias_point(point_s32 p) {
    return {bias_coord(p.x()), bias_coord(p.y())};
}
inline point_s32 unbias_point(point p) {
    return {unbias_coord(p.x()), unbias_coord(p.y())};
}

inline ring bias_ring(const ring_s32& src) {
    ring r;
    for (const auto& pt : src) r.push_back(bias_point(pt));
    return r;
}
inline ring_s32 unbias_ring(const ring& src) {
    ring_s32 r;
    for (const auto& pt : src) r.push_back(unbias_point(pt));
    return r;
}

inline polygon bias_poly(const polygon_s32& src) {
    polygon p;
    p.outer() = bias_ring(src.outer());
    for (const auto& inner : src.inners()) p.inners().push_back(bias_ring(inner));
    return p;
}
inline polygon_s32 unbias_poly(const polygon& src) {
    polygon_s32 p;
    p.outer() = unbias_ring(src.outer());
    for (const auto& inner : src.inners()) p.inners().push_back(unbias_ring(inner));
    return p;
}

inline multi_polygon_u32 bias_mp(const multi_polygon_s32& src) {
    multi_polygon_u32 mp;
    for (const auto& poly : src) mp.push_back(bias_poly(poly));
    return mp;
}
inline multi_polygon_s32 unbias_mp(const multi_polygon_u32& src) {
    multi_polygon_s32 mp;
    for (const auto& poly : src) mp.push_back(unbias_poly(poly));
    return mp;
}

inline multi_polygon_s32 add(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    return unbias_mp(best_clipper::add(bias_mp(a), bias_mp(b)));
}

inline multi_polygon_s32 intersection(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    return unbias_mp(best_clipper::intersection(bias_mp(a), bias_mp(b)));
}

inline multi_polygon_s32 xor_(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    return unbias_mp(best_clipper::xor_(bias_mp(a), bias_mp(b)));
}

inline multi_polygon_s32 difference(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    return unbias_mp(best_clipper::difference(bias_mp(a), bias_mp(b)));
}

inline multi_polygon_s32 self_or(const multi_polygon_s32& a) {
    return unbias_mp(best_clipper::self_or(bias_mp(a)));
}
