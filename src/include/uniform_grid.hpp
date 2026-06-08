#pragma once

/// Uniform grid spatial index for int32_t-coordinate 2D bounding boxes.
///
/// Divides space into fixed-size cells. Each item (box + value) is placed
/// into every cell its box touches. Uses sparse storage — only non-empty
/// cells consume memory. Falls back to linear scan for fewer than THRESHOLD
/// items.

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include <boost/geometry.hpp>

namespace bg = boost::geometry;

namespace best_clipper::uniform_grid {

namespace detail {
constexpr size_t THRESHOLD = 256;

inline int64_t to_int(int32_t v) { return (int64_t)v; }

inline int64_t cell_key(int64_t cx, int64_t cy, int64_t x_cells) {
    return cy * x_cells + cx;
}
}  // namespace detail

template <typename Value>
struct grid {
    using box = bg::model::box<bg::model::d2::point_xy<int32_t>>;

    grid(int64_t cell_size = 2) : _min_x(0), _min_y(0), _x_cells(0), _y_cells(0), _cell_size(cell_size) {}

    /// Build from [ (box, value), ... ]. Pass cell_size=0 to auto-compute.
    template <typename Range>
    explicit grid(Range&& items, int64_t cell_size = 2) : _cell_size(cell_size) {
        auto n = items.size();
        _values.reserve(n);
        _boxes.reserve(n);
        for (auto&& [b, v] : items) {
            _boxes.push_back(b);
            _values.push_back(std::move(v));
        }

        if (n < detail::THRESHOLD) {
            _x_cells = 0;  // signal linear-scan mode
            _y_cells = 0;
            return;
        }

        // Compute bounds
        int64_t min_x = INT64_MAX, min_y = INT64_MAX, max_x = 0, max_y = 0;
        for (auto& b : _boxes) {
            auto x0 = detail::to_int(bg::get<0, 0>(b));
            auto y0 = detail::to_int(bg::get<0, 1>(b));
            auto x1 = detail::to_int(bg::get<1, 0>(b));
            auto y1 = detail::to_int(bg::get<1, 1>(b));
            if (x0 < min_x) min_x = x0; if (y0 < min_y) min_y = y0;
            if (x1 > max_x) max_x = x1; if (y1 > max_y) max_y = y1;
        }
        _min_x = min_x;
        _min_y = min_y;
        if (_cell_size == 0) {
            int64_t dx = max_x - min_x;
            int64_t dy = max_y - min_y;
            int64_t cells_per_dim = (int64_t)std::sqrt((double)n);
            if (cells_per_dim < 1) cells_per_dim = 1;
            _cell_size = std::max((int64_t)1, std::max(dx, dy) / cells_per_dim);
        }
        _x_cells = (max_x - min_x) / _cell_size + 1;
        _y_cells = (max_y - min_y) / _cell_size + 1;

        // Assign items to cells: collect (cell_key, item_index) pairs, sort, group
        std::vector<std::pair<int64_t, size_t>> cell_pairs;
        cell_pairs.reserve(_boxes.size() * 4);  // rough estimate
        for (size_t i = 0; i < _boxes.size(); i++) {
            auto& b = _boxes[i];
            int64_t cx1 = (detail::to_int(bg::get<0, 0>(b)) - min_x) / _cell_size;
            int64_t cy1 = (detail::to_int(bg::get<0, 1>(b)) - min_y) / _cell_size;
            int64_t cx2 = (detail::to_int(bg::get<1, 0>(b)) - min_x) / _cell_size;
            int64_t cy2 = (detail::to_int(bg::get<1, 1>(b)) - min_y) / _cell_size;
            if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
            if (cx2 >= _x_cells) cx2 = _x_cells - 1;
            if (cy2 >= _y_cells) cy2 = _y_cells - 1;
            for (int64_t cy = cy1; cy <= cy2; cy++)
                for (int64_t cx = cx1; cx <= cx2; cx++)
                    cell_pairs.emplace_back(detail::cell_key(cx, cy, _x_cells), i);
        }

        // Sort by cell key, then group into contiguous arrays
        std::sort(cell_pairs.begin(), cell_pairs.end(),
            [](auto& a, auto& b) { return a.first < b.first; });

        _cell_keys.reserve(cell_pairs.size() / 2);  // rough
        _cell_begins.reserve(cell_pairs.size() / 2 + 1);
        _cell_begins.push_back(0);
        _cell_items.reserve(cell_pairs.size());
        _cell_max_x.reserve(cell_pairs.size() / 2);
        for (size_t i = 0; i < cell_pairs.size(); ) {
            auto key = cell_pairs[i].first;
            _cell_keys.push_back(key);
            size_t j = i;
            int64_t max_x = 0;
            while (j < cell_pairs.size() && cell_pairs[j].first == key) {
                size_t item_idx = cell_pairs[j].second;
                _cell_items.push_back(item_idx);
                auto x1 = detail::to_int(bg::get<0, 0>(_boxes[item_idx]));
                auto x2 = detail::to_int(bg::get<1, 0>(_boxes[item_idx]));
                if (x2 > max_x) max_x = x2;
                j++;
            }
            _cell_max_x.push_back(max_x);
            _cell_begins.push_back(_cell_items.size());
            i = j;
        }
    }

    size_t size() const { return _values.size(); }

    /// Calls cb(value_index) for each value whose box intersects query_box.
    template <typename Callback>
    void query_intersects(const box& query_box, Callback&& cb) const {
        query_impl(query_box, std::forward<Callback>(cb));
    }

private:
    std::vector<box> _boxes;
    std::vector<Value> _values;
    int64_t _min_x, _min_y, _x_cells, _y_cells, _cell_size;
    std::vector<int64_t> _cell_keys;     // sorted cell key per non-empty cell
    std::vector<size_t> _cell_begins;    // begin offset per cell (+1 sentinel)
    std::vector<size_t> _cell_items;     // flat item indices
    std::vector<int64_t> _cell_max_x;    // max box max_x per cell (for fast skip)

    template <typename Callback>
    void query_impl(const box& query_box, Callback&& cb) const {
        int64_t qx0 = detail::to_int(bg::get<0, 0>(query_box));
        int64_t qy0 = detail::to_int(bg::get<0, 1>(query_box));
        int64_t qx1 = detail::to_int(bg::get<1, 0>(query_box));
        int64_t qy1 = detail::to_int(bg::get<1, 1>(query_box));

        if (_x_cells == 0) {
            // Linear-scan mode
            for (size_t i = 0; i < _boxes.size(); i++) {
                auto& b = _boxes[i];
                if (qx1 < detail::to_int(bg::get<0, 0>(b))
                    || qx0 > detail::to_int(bg::get<1, 0>(b))
                    || qy1 < detail::to_int(bg::get<0, 1>(b))
                    || qy0 > detail::to_int(bg::get<1, 1>(b)))
                    continue;
                cb(i);
            }
            return;
        }

        // Grid mode
        int64_t cx1 = (qx0 - _min_x) / _cell_size;
        int64_t cy1 = (qy0 - _min_y) / _cell_size;
        int64_t cx2 = (qx1 - _min_x) / _cell_size;
        int64_t cy2 = (qy1 - _min_y) / _cell_size;
        if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
        if (cx2 >= _x_cells) cx2 = _x_cells - 1;
        if (cy2 >= _y_cells) cy2 = _y_cells - 1;

        // Single-row optimization: keys are consecutive, iterate right-to-left
        // so callers can use a max_x filter to skip most cells early.
        if (cy1 == cy2) {
            int64_t base = detail::cell_key(0, cy1, _x_cells);
            auto it = std::lower_bound(_cell_keys.begin(), _cell_keys.end(), base + cx1);
            auto end_it = std::lower_bound(it, _cell_keys.end(), base + cx2 + 1);
            // Iterate backwards: rightmost cells first
            while (end_it != it) {
                --end_it;
                size_t idx = (size_t)(end_it - _cell_keys.begin());
                for (auto j = _cell_begins[idx]; j != _cell_begins[idx + 1]; ++j) {
                    size_t i = _cell_items[j];
                    auto& b = _boxes[i];
                    if (qx1 < detail::to_int(bg::get<0, 0>(b))
                        || qx0 > detail::to_int(bg::get<1, 0>(b))
                        || qy1 < detail::to_int(bg::get<0, 1>(b))
                        || qy0 > detail::to_int(bg::get<1, 1>(b)))
                        continue;
                    cb(i);
                }
            }
            return;
        }

        for (int64_t cy = cy1; cy <= cy2; cy++) {
            int64_t base = detail::cell_key(0, cy, _x_cells);
            auto it = std::lower_bound(_cell_keys.begin(), _cell_keys.end(), base + cx1);
            auto end_it = std::lower_bound(it, _cell_keys.end(), base + cx2 + 1);
            for (; it != end_it; ++it) {
                size_t idx = (size_t)(it - _cell_keys.begin());
                for (auto j = _cell_begins[idx]; j != _cell_begins[idx + 1]; ++j) {
                    size_t i = _cell_items[j];
                    auto& b = _boxes[i];
                    if (qx1 < detail::to_int(bg::get<0, 0>(b))
                        || qx0 > detail::to_int(bg::get<1, 0>(b))
                        || qy1 < detail::to_int(bg::get<0, 1>(b))
                        || qy0 > detail::to_int(bg::get<1, 1>(b)))
                        continue;
                    cb(i);
                }
            }
        }
    }
};

}  // namespace best_clipper::uniform_grid
