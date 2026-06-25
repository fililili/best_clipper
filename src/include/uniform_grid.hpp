#pragma once

/// Uniform grid spatial index for coordinate_type 2D bounding boxes.
///
/// Divides space into fixed-size cells. Each item (box + index) is placed
/// into every cell its box touches. The grid only stores cell→item-index
/// mappings; callers retain boxes and pass them to queries for filtering.

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "geometry_types.hpp"

namespace best_clipper::uniform_grid {

struct grid {

  std::size_t get_flat_index(coordinate_type x_cell,
                             coordinate_type y_cell) const {
    assert(x_cell >= 0 && x_cell < _x_cells);
    assert(y_cell >= 0 && y_cell < _y_cells);
    return (std::size_t)y_cell * _x_cells + x_cell;
  }

  coordinate_type cell_x(coordinate_type x) const {
    return std::clamp(((x - _min_x) / _cell_size), (coordinate_type)0,
                      _x_cells - 1);
  }
  coordinate_type cell_y(coordinate_type y) const {
    return std::clamp(((y - _min_y) / _cell_size), (coordinate_type)0,
                      _y_cells - 1);
  }

  std::size_t cells() { return (std::size_t)_x_cells * _y_cells; }

  /// Build from boxes. Item index = position in the vector.
  explicit grid(const std::vector<box> &items) {
    auto n = items.size();

    coordinate_type min_x = std::numeric_limits<coordinate_type>::max();
    coordinate_type min_y = std::numeric_limits<coordinate_type>::max();
    coordinate_type max_x = std::numeric_limits<coordinate_type>::min();
    coordinate_type max_y = std::numeric_limits<coordinate_type>::min();
    for (auto &b : items) {
      auto x0 = bg::get<0, 0>(b), y0 = bg::get<0, 1>(b);
      auto x1 = bg::get<1, 0>(b), y1 = bg::get<1, 1>(b);
      if (x0 < min_x)
        min_x = x0;
      if (y0 < min_y)
        min_y = y0;
      if (x1 > max_x)
        max_x = x1;
      if (y1 > max_y)
        max_y = y1;
    }
    _min_x = min_x;
    _min_y = min_y;
    coordinate_type dx = max_x - min_x;
    coordinate_type dy = max_y - min_y;
    auto cells_per_dim = (coordinate_type)std::sqrt((double)n);
    if (cells_per_dim < 1)
      cells_per_dim = 1;
    _cell_size = std::max((coordinate_type)1, std::max(dx, dy) / cells_per_dim);
    _x_cells = (max_x - min_x) / _cell_size + 1;
    _y_cells = (max_y - min_y) / _cell_size + 1;
    std::size_t num_cells = (std::size_t)_x_cells * _y_cells;

    std::vector<std::pair<std::size_t, std::size_t>> cell_pairs;
    cell_pairs.reserve(n * 4);
    for (std::size_t i = 0; i < n; i++) {
      auto &b = items[i];
      coordinate_type bx0 = bg::get<0, 0>(b), bx1 = bg::get<1, 0>(b);
      coordinate_type by0 = bg::get<0, 1>(b), by1 = bg::get<1, 1>(b);
      if (bx0 > bx1 || by0 > by1)
        continue;
      coordinate_type cx1 = cell_x(bx0);
      coordinate_type cy1 = cell_y(by0);
      coordinate_type cx2 = cell_x(bx1);
      coordinate_type cy2 = cell_y(by1);
      for (coordinate_type cy = cy1; cy <= cy2; cy++)
        for (coordinate_type cx = cx1; cx <= cx2; cx++)
          cell_pairs.emplace_back(get_flat_index(cx, cy), i);
    }

    std::vector<size_t> cell_counts(num_cells, 0);
    for (auto &[cell, item] : cell_pairs)
      cell_counts[cell]++;

    std::vector<size_t> begins(num_cells + 1, 0);
    for (size_t c = 0; c < num_cells; c++)
      begins[c + 1] = begins[c] + cell_counts[c];

    _cell_items.resize(cell_pairs.size());
    auto cursors = begins;
    for (auto &[cell, item] : cell_pairs)
      _cell_items[cursors[cell]++] = item;

    _cell_begins = std::move(begins);
  }

  template <typename Callback>
  void query_intersects(const box &query_box, Callback &&cb) const {
    coordinate_type qx0 = bg::get<0, 0>(query_box);
    coordinate_type qy0 = bg::get<0, 1>(query_box);
    coordinate_type qx1 = bg::get<1, 0>(query_box);
    coordinate_type qy1 = bg::get<1, 1>(query_box);

    coordinate_type cx1 = cell_x(qx0);
    coordinate_type cy1 = cell_y(qy0);
    coordinate_type cx2 = cell_x(qx1);
    coordinate_type cy2 = cell_y(qy1);

    for (coordinate_type cy = cy1; cy <= cy2; cy++) {
      for (coordinate_type cx = cx1; cx <= cx2; cx++) {
        size_t ci = get_flat_index(cx, cy);
        for (auto j = _cell_begins[ci]; j != _cell_begins[ci + 1]; ++j)
          cb(_cell_items[j]);
      }
    }
  }

  template <typename Callback>
  void query_ray_left(coordinate_type query_x, coordinate_type query_y,
                      Callback &&cb) const {
    coordinate_type cx = cell_x(query_x);
    coordinate_type cy = cell_y(query_y);
    coordinate_type best_x_floor = std::numeric_limits<coordinate_type>::min();
    // best_x_real is not interger, so we store the floor of it. best_x_real >=
    // best_x_floor
    for (coordinate_type current_cx = cx; current_cx != (coordinate_type)(-1);
         current_cx--) {
      if (best_x_floor >
          _min_x + (multi_coordinate_type)(current_cx + 1) * _cell_size) {
        break;
      }
      std::size_t ci = get_flat_index(current_cx, cy);
      for (auto j = _cell_begins[ci]; j != _cell_begins[ci + 1]; ++j) {
        best_x_floor = cb(_cell_items[j]);
      }
    }
  }

  coordinate_type _min_x = 0, _min_y = 0, _cell_size = 0;
  coordinate_type _x_cells = 0; // use coordinate_type because _x_cells < max_x
  coordinate_type _y_cells = 0;
  std::vector<size_t>
      _cell_begins; // dense: (x_cells * y_cells + 1) offsets into _cell_items
  std::vector<size_t> _cell_items; // flat item indices, grouped by cell
};

} // namespace best_clipper::uniform_grid
