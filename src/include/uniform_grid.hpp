#pragma once

/// Uniform grid spatial index for int32_t-coordinate 2D bounding boxes.
///
/// Divides space into fixed-size cells. Each item (box + index) is placed
/// into every cell its box touches. The grid only stores cell→item-index
/// mappings; callers retain boxes and pass them to queries for filtering.

#include <cstdint>
#include <utility>
#include <vector>

#include <boost/geometry.hpp>

namespace bg = boost::geometry;

namespace best_clipper::uniform_grid {

struct grid {
    using box = bg::model::box<bg::model::d2::point_xy<int32_t>>;

    grid() : _min_x(0), _min_y(0), _x_cells(0), _y_cells(0), _cell_size(0) {}
    explicit grid(int32_t cell_size) : _min_x(0), _min_y(0), _x_cells(0), _y_cells(0), _cell_size(cell_size) {}

    /// Build from [ (box, index), ... ]. Pass cell_size=0 to auto-compute.
    explicit grid(const std::vector<std::pair<box, std::size_t>>& items, int32_t cell_size = 0)
        : _cell_size(cell_size) {
        auto n = items.size();
        if (n == 0) return;

        // Compute bounds (skip invalid boxes)
        int32_t min_x = INT32_MAX, min_y = INT32_MAX, max_x = INT32_MIN, max_y = INT32_MIN;
        for (auto& [b, v] : items) {
            auto x0 = bg::get<0, 0>(b), y0 = bg::get<0, 1>(b);
            auto x1 = bg::get<1, 0>(b), y1 = bg::get<1, 1>(b);
            if (x0 > x1 || y0 > y1) continue;
            if (x0 < min_x) min_x = x0; if (y0 < min_y) min_y = y0;
            if (x1 > max_x) max_x = x1; if (y1 > max_y) max_y = y1;
        }
        _min_x = min_x;
        _min_y = min_y;
        if (_cell_size == 0) {
            int32_t dx = max_x - min_x;
            int32_t dy = max_y - min_y;
            int32_t cells_per_dim = (int32_t)std::sqrt((double)n);
            if (cells_per_dim < 1) cells_per_dim = 1;
            _cell_size = std::max((int32_t)1, std::max(dx, dy) / cells_per_dim);
        }
        _x_cells = (max_x - min_x) / _cell_size + 1;
        _y_cells = (max_y - min_y) / _cell_size + 1;
        size_t num_cells = (size_t)(_x_cells * _y_cells);

        // Collect (cell_index, item_index) pairs
        std::vector<std::pair<size_t, size_t>> cell_pairs;
        cell_pairs.reserve(n * 4);
        for (size_t i = 0; i < n; i++) {
            auto& b = items[i].first;
            int32_t bx0 = bg::get<0, 0>(b), bx1 = bg::get<1, 0>(b);
            int32_t by0 = bg::get<0, 1>(b), by1 = bg::get<1, 1>(b);
            if (bx0 > bx1 || by0 > by1) continue;
            int32_t cx1 = (bx0 - min_x) / _cell_size;
            int32_t cy1 = (by0 - min_y) / _cell_size;
            int32_t cx2 = (bx1 - min_x) / _cell_size;
            int32_t cy2 = (by1 - min_y) / _cell_size;
            if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
            if (cx2 >= _x_cells) cx2 = _x_cells - 1;
            if (cy2 >= _y_cells) cy2 = _y_cells - 1;
            for (int32_t cy = cy1; cy <= cy2; cy++)
                for (int32_t cx = cx1; cx <= cx2; cx++)
                    cell_pairs.emplace_back((size_t)(cy * _x_cells + cx), i);
        }

        // Counting sort by cell index
        std::vector<size_t> cell_counts(num_cells, 0);
        for (auto& [cell, item] : cell_pairs)
            cell_counts[cell]++;

        std::vector<size_t> begins(num_cells + 1, 0);
        for (size_t c = 0; c < num_cells; c++)
            begins[c + 1] = begins[c] + cell_counts[c];

        _cell_items.resize(cell_pairs.size());
        auto cursors = begins;
        for (auto& [cell, item] : cell_pairs) {
            _cell_items[cursors[cell]++] = item;
        }

        _cell_begins = std::move(begins);
    }

    template <typename Callback>
    void query_intersects(const box& query_box, Callback&& cb) const {
        int32_t qx0 = bg::get<0, 0>(query_box);
        int32_t qy0 = bg::get<0, 1>(query_box);
        int32_t qx1 = bg::get<1, 0>(query_box);
        int32_t qy1 = bg::get<1, 1>(query_box);

        if (_cell_begins.empty()) return;

        int32_t cx1 = (qx0 - _min_x) / _cell_size;
        int32_t cy1 = (qy0 - _min_y) / _cell_size;
        int32_t cx2 = (qx1 - _min_x) / _cell_size;
        int32_t cy2 = (qy1 - _min_y) / _cell_size;
        if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
        if (cx2 >= _x_cells) cx2 = _x_cells - 1;
        if (cy2 >= _y_cells) cy2 = _y_cells - 1;

        if (cy1 == cy2) {
            for (int32_t cx = cx2; cx >= cx1; cx--) {
                size_t ci = (size_t)(cy1 * _x_cells + cx);
                for (auto j = _cell_begins[ci]; j != _cell_begins[ci + 1]; ++j)
                    cb(_cell_items[j]);
            }
            return;
        }

        for (int32_t cy = cy1; cy <= cy2; cy++) {
            for (int32_t cx = cx1; cx <= cx2; cx++) {
                size_t ci = (size_t)(cy * _x_cells + cx);
                for (auto j = _cell_begins[ci]; j != _cell_begins[ci + 1]; ++j)
                    cb(_cell_items[j]);
            }
        }
    }

    int32_t _min_x = 0, _min_y = 0, _x_cells = 0, _y_cells = 0, _cell_size = 0;
    std::vector<size_t> _cell_begins;  // dense: (x_cells * y_cells + 1) offsets into _cell_items
    std::vector<size_t> _cell_items;   // flat item indices, grouped by cell
};

}  // namespace best_clipper::uniform_grid
