#pragma once

#include <boost/geometry.hpp>

namespace best_clipper {

namespace bg = boost::geometry;
using coordinate_type = int32_t;
using point = bg::model::d2::point_xy<coordinate_type>;
using segment = bg::model::segment<point>;
using box = bg::model::box<point>;
using ring = bg::model::ring<point>;
using polygon = bg::model::polygon<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;
} // namespace best_clipper