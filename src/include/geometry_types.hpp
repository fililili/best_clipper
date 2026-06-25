#pragma once

#include <boost/geometry.hpp>
#include <cstdint>

namespace best_clipper {

namespace bg = boost::geometry;
// Addition/subtraction do not guard against overflow.
// Multiplication results use multi_coordinate_type.
using coordinate_type = std::int32_t;
using multi_coordinate_type = std::int64_t;
using point = bg::model::d2::point_xy<coordinate_type>;
using segment = bg::model::segment<point>;
using box = bg::model::box<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;

inline bool bbox_overlap(const box &a, const box &b) {
  return bg::get<0, 0>(a) <= bg::get<1, 0>(b) &&
         bg::get<1, 0>(a) >= bg::get<0, 0>(b) &&
         bg::get<0, 1>(a) <= bg::get<1, 1>(b) &&
         bg::get<1, 1>(a) >= bg::get<0, 1>(b);
}

} // namespace best_clipper