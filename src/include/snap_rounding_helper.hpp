#pragma once

#include "geometry_types.hpp"
#include <cstdint>
#include <optional>
#include <utility>

namespace best_clipper {

// Order hot pixels along a segment.  Every pixel passed to operator() lies on
// the segment, so p = (x1 + t*dx, y1 + t*dy) for some t in [0, 1].  The dot
// product with (dx, dy) is x1*dx + y1*dy + t*(dx*dx+dy*dy) - monotonic in t.
// Replacing (dx, dy) by their signs (sx, sy) in {-1,0,1} gives:
//   dot(p, (sx, sy)) = x1*sx + y1*sy + t*(|dx| + |dy|)
// The constant offset cancels in comparison and |dx|+|dy| > 0, so the ordering
// by t is preserved while the comparison simplifies to sign-flip arithmetic.
struct less_by_segment {
  int sx, sy;

  less_by_segment(const segment &s) {
    coordinate_type d = bg::get<1, 0>(s) - bg::get<0, 0>(s);
    sx = (d > 0) - (d < 0);
    d = bg::get<1, 1>(s) - bg::get<0, 1>(s);
    sy = (d > 0) - (d < 0);
  }

  bool operator()(point p1, point p2) const {
    return bg::get<0>(p1) * sx + bg::get<1>(p1) * sy <
           bg::get<0>(p2) * sx + bg::get<1>(p2) * sy;
  }
};

inline multi_coordinate_type cross(coordinate_type dx1, coordinate_type dy1,
                                   coordinate_type dx2, coordinate_type dy2) {
  return (multi_coordinate_type)dx1 * dy2 - (multi_coordinate_type)dy1 * dx2;
}

// Hot pixel Pi(p) = (x-0.5, x+0.5] x (y-0.5, y+0.5] - unit square centered at
// integral point p.  Right/top closed, left/bottom open.  This half-open
// convention guarantees each plane point belongs to exactly one pixel and is
// the foundation of the topological consistency proof in Guibas & Marimont
// (1998, "Rounding Arrangements Dynamically", section 2-3).
// Snap rounding in get_intersection() rounds half-integers toward -inf
// (left/bottom), consistent with the open boundaries.
inline bool is_point_on_segment(point p, segment s) {
  coordinate_type x = bg::get<0>(p), y = bg::get<1>(p);
  coordinate_type x1 = bg::get<0, 0>(s), y1 = bg::get<0, 1>(s);
  coordinate_type x2 = bg::get<1, 0>(s), y2 = bg::get<1, 1>(s);
  coordinate_type dx = x2 - x1, dy = y2 - y1;
  coordinate_type adx = dx >= 0 ? dx : -dx;
  coordinate_type ady = dy >= 0 ? dy : -dy;
  coordinate_type T = adx + ady;
  if (T == 0)
    return x == x1 && y == y1;

  {
    multi_coordinate_type dot = (multi_coordinate_type)(x - x1) * dx +
                                (multi_coordinate_type)(y - y1) * dy;
    multi_coordinate_type lhs = dot + dot + T;
    if (dx >= 0 && dy >= 0) {
      if (lhs < 0)
        return false;
    } else {
      if (lhs <= 0)
        return false;
    }
  }
  {
    multi_coordinate_type dot = (multi_coordinate_type)(x2 - x) * dx +
                                (multi_coordinate_type)(y2 - y) * dy;
    multi_coordinate_type lhs = dot + dot + T;
    if (dx <= 0 && dy <= 0) {
      if (lhs < 0)
        return false;
    } else {
      if (lhs <= 0)
        return false;
    }
  }
  {
    multi_coordinate_type cross = (multi_coordinate_type)(x - x1) * dy -
                                  (multi_coordinate_type)(y - y1) * dx;
    multi_coordinate_type abs2 =
        cross >= 0 ? cross + cross : (-cross) + (-cross);
    if (cross >= 0) {
      if (dx >= 0 && dy <= 0) {
        if (abs2 > T)
          return false;
      } else {
        if (abs2 >= T)
          return false;
      }
    } else {
      if (dx <= 0 && dy >= 0) {
        if (abs2 > T)
          return false;
      } else {
        if (abs2 >= T)
          return false;
      }
    }
  }
  return true;
}

inline std::optional<point> get_intersection(segment s1, segment s2) {
  coordinate_type x1 = bg::get<0, 0>(s1), y1 = bg::get<0, 1>(s1);
  coordinate_type x2 = bg::get<1, 0>(s1), y2 = bg::get<1, 1>(s1);
  coordinate_type x3 = bg::get<0, 0>(s2), y3 = bg::get<0, 1>(s2);
  coordinate_type x4 = bg::get<1, 0>(s2), y4 = bg::get<1, 1>(s2);
  coordinate_type dx1 = x2 - x1, dy1 = y2 - y1;
  coordinate_type dx2 = x4 - x3, dy2 = y4 - y3;

  // Fast path: axis-aligned segments.
  if (dx1 == 0 && dy2 == 0) {
    if (x1 <= std::min(x3, x4) || x1 >= std::max(x3, x4))
      return {};
    if (y3 <= std::min(y1, y2) || y3 >= std::max(y1, y2))
      return {};
    return point{x1, y3};
  }
  if (dy1 == 0 && dx2 == 0) {
    if (x3 <= std::min(x1, x2) || x3 >= std::max(x1, x2))
      return {};
    if (y1 <= std::min(y3, y4) || y1 >= std::max(y3, y4))
      return {};
    return point{x3, y1};
  }

  if (cross(dx1, dy1, dx2, dy2) == 0)
    return {};

  multi_coordinate_type d = cross(dx1, dy1, dx2, dy2);
  multi_coordinate_type t1_num = cross(x3 - x1, y3 - y1, dx2, dy2);
  multi_coordinate_type t2_num = cross(x3 - x1, y3 - y1, dx1, dy1);
  if (d < 0) {
    d = -d;
    t1_num = -t1_num;
    t2_num = -t2_num;
  }
  if (t1_num <= 0 || t1_num >= d || t2_num <= 0 || t2_num >= d)
    return {};

  auto snap = [&](coordinate_type x, coordinate_type dx,
                  multi_coordinate_type num,
                  multi_coordinate_type den) -> coordinate_type {
    multi_coordinate_type num_dx = num * dx;
    multi_coordinate_type n = 2 * num_dx + den - 1;
    multi_coordinate_type dd = 2 * den;
    multi_coordinate_type q = n >= 0 ? n / dd : (n - dd + 1) / dd;
    return (coordinate_type)(x + q);
  };
  auto ret = point{snap(x1, dx1, t1_num, d), snap(y1, dy1, t1_num, d)};
  assert(is_point_on_segment(ret, s1));
  assert(is_point_on_segment(ret, s2));
  return ret;
}

} // namespace best_clipper
