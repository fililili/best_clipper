#pragma once
#include "geometry_types.hpp"

namespace best_clipper {

multi_polygon union_(const multi_polygon &a, const multi_polygon &b);
multi_polygon intersection(const multi_polygon &a, const multi_polygon &b);
multi_polygon symmetric_difference(const multi_polygon &a,
                                   const multi_polygon &b);
multi_polygon difference(const multi_polygon &a, const multi_polygon &b);
multi_polygon robust_self_or(const multi_polygon &a);

} // namespace best_clipper
