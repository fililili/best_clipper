#pragma once

/// Uniform grid spatial index for int32_t-coordinate 2D bounding boxes.
///
/// Divides space into fixed-size cells. Each item (box + value) is placed
/// into every cell its box touches.

#include <cstdint>
#include <utility>
#include <vector>

#include <boost/geometry.hpp>

namespace bg = boost::geometry;

namespace best_clipper::uniform_grid {

template <typename Value>
struct grid {
    using box = bg::model::box<bg::model::d2::point_xy<int32_t>>;

    grid(int32_t cell_size) : _min_x(0), _min_y(0), _x_cells(0), _y_cells(0), _cell_size(cell_size) {}

    /// Build from [ (box, value), ... ]. Pass cell_size=0 to auto-compute.
    template <typename Range>
    explicit grid(Range&& items, int32_t cell_size) : _cell_size(cell_size) {
        auto n = items.size();
        _values.reserve(n);
        _boxes.reserve(n);
        for (auto&& [b, v] : items) {
            _boxes.push_back(b);
            _values.push_back(std::move(v));
        }

        if (n == 0) return;

        // Compute bounds
        int32_t min_x = INT32_MAX, min_y = INT32_MAX, max_x = INT32_MIN, max_y = INT32_MIN;
        for (auto& b : _boxes) {
            auto x0 = bg::get<0, 0>(b);
            auto y0 = bg::get<0, 1>(b);
            auto x1 = bg::get<1, 0>(b);
            auto y1 = bg::get<1, 1>(b);
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

        // Assign items to cells: collect (cell_index, item_index) pairs, sort, group
        std::vector<std::pair<size_t, size_t>> cell_pairs;
        cell_pairs.reserve(_boxes.size() * 4);
        for (size_t i = 0; i < _boxes.size(); i++) {
            auto& b = _boxes[i];
            int32_t cx1 = (bg::get<0, 0>(b) - min_x) / _cell_size;
            int32_t cy1 = (bg::get<0, 1>(b) - min_y) / _cell_size;
            int32_t cx2 = (bg::get<1, 0>(b) - min_x) / _cell_size;
            int32_t cy2 = (bg::get<1, 1>(b) - min_y) / _cell_size;
            if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
            if (cx2 >= _x_cells) cx2 = _x_cells - 1;
            if (cy2 >= _y_cells) cy2 = _y_cells - 1;
            for (int32_t cy = cy1; cy <= cy2; cy++)
                for (int32_t cx = cx1; cx <= cx2; cx++)
                    cell_pairs.emplace_back((size_t)(cy * _x_cells + cx), i);
        }

        // Assign items to cells using counting sort (bucket sort by cell index)
        std::vector<size_t> cell_counts(num_cells, 0);
        struct item_cell { size_t cell; size_t item; };
        std::vector<item_cell> temp(cell_pairs.size());
        for (size_t i = 0; i < cell_pairs.size(); i++) {
            temp[i] = {cell_pairs[i].first, cell_pairs[i].second};
            cell_counts[cell_pairs[i].first]++;
        }

        // Free cell_pairs early — no longer needed
        cell_pairs.clear();
        cell_pairs.shrink_to_fit();

        std::vector<size_t> cell_begins_temp(num_cells + 1, 0);
        for (size_t c = 0; c < num_cells; c++)
            cell_begins_temp[c + 1] = cell_begins_temp[c] + cell_counts[c];

        _cell_items.resize(temp.size());
        _cell_max_x.assign(num_cells, INT32_MIN);
        auto cursors = cell_begins_temp;
        for (auto& ic : temp) {
            auto pos = cursors[ic.cell]++;
            _cell_items[pos] = ic.item;
            auto x2 = bg::get<1, 0>(_boxes[ic.item]);
            if (x2 > _cell_max_x[ic.cell]) _cell_max_x[ic.cell] = x2;
        }
        temp.clear();

        _cell_begins = std::move(cell_begins_temp);
    }

    size_t size() const { return _values.size(); }

    template <typename Callback>
    void query_intersects(const box& query_box, Callback&& cb) const {
        int32_t qx0 = bg::get<0, 0>(query_box);
        int32_t qy0 = bg::get<0, 1>(query_box);
        int32_t qx1 = bg::get<1, 0>(query_box);
        int32_t qy1 = bg::get<1, 1>(query_box);

        if (_boxes.empty()) return;

        int32_t cx1 = (qx0 - _min_x) / _cell_size;
        int32_t cy1 = (qy0 - _min_y) / _cell_size;
        int32_t cx2 = (qx1 - _min_x) / _cell_size;
        int32_t cy2 = (qy1 - _min_y) / _cell_size;
        if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
        if (cx2 >= _x_cells) cx2 = _x_cells - 1;
        if (cy2 >= _y_cells) cy2 = _y_cells - 1;

        // Single-row: iterate backwards so callers can skip early via max_x.
        if (cy1 == cy2) {
            for (int32_t cx = cx2; cx >= cx1; cx--) {
                size_t ci = (size_t)(cy1 * _x_cells + cx);
                for (auto j = _cell_begins[ci]; j != _cell_begins[ci + 1]; ++j) {
                    size_t i = _cell_items[j];
                    auto& b = _boxes[i];
                    if (qx1 < bg::get<0, 0>(b)
                        || qx0 > bg::get<1, 0>(b)
                        || qy1 < bg::get<0, 1>(b)
                        || qy0 > bg::get<1, 1>(b))
                        continue;
                    cb(i);
                }
            }
            return;
        }

        for (int32_t cy = cy1; cy <= cy2; cy++) {
            for (int32_t cx = cx1; cx <= cx2; cx++) {
                size_t ci = (size_t)(cy * _x_cells + cx);
                for (auto j = _cell_begins[ci]; j != _cell_begins[ci + 1]; ++j) {
                    size_t i = _cell_items[j];
                    auto& b = _boxes[i];
                    if (qx1 < bg::get<0, 0>(b)
                        || qx0 > bg::get<1, 0>(b)
                        || qy1 < bg::get<0, 1>(b)
                        || qy0 > bg::get<1, 1>(b))
                        continue;
                    cb(i);
                }
            }
        }
    }

private:
    std::vector<box> _boxes;
    std::vector<Value> _values;
    int32_t _min_x, _min_y, _x_cells, _y_cells, _cell_size;
    std::vector<size_t> _cell_begins;  // dense: (x_cells * y_cells + 1) offsets into _cell_items
    std::vector<size_t> _cell_items;   // flat item indices, grouped by cell
    std::vector<int32_t> _cell_max_x;  // max box max_x per cell
};

}  // namespace best_clipper::uniform_grid
