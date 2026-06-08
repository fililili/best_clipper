#include <gtest/gtest.h>
#include <boost/geometry.hpp>
#include "uniform_grid.hpp"

namespace bg = boost::geometry;
using point = bg::model::d2::point_xy<int32_t>;
using box = bg::model::box<point>;
using best_clipper::uniform_grid::grid;

// Helper: build grid from vector of boxes
template <typename BoxVec>
static auto build_grid(BoxVec&& boxes) {
    std::vector<std::pair<box, size_t>> items(boxes.size());
    for (size_t i = 0; i < boxes.size(); i++)
        items[i] = {boxes[i], i};
    return grid<size_t>(std::move(items), 0);
}

TEST(UniformGrid, Empty) {
    grid<size_t> g(2);
    EXPECT_EQ(g.size(), 0u);
    box q{point{0, 0}, point{100, 100}};
    int count = 0;
    g.query_intersects(q, [&](size_t) { count++; });
    EXPECT_EQ(count, 0);
}

TEST(UniformGrid, SingleItem) {
    std::vector<box> boxes = {{point{0, 0}, point{10, 10}}};
    auto g = build_grid(boxes);

    // Point inside
    int count = 0;
    g.query_intersects(box{point{5, 5}, point{5, 5}}, [&](size_t i) {
        EXPECT_EQ(i, 0u);
        count++;
    });
    EXPECT_EQ(count, 1);
}

TEST(UniformGrid, TwoDisjoint) {
    std::vector<box> boxes = {
        {point{0, 0}, point{10, 10}},
        {point{100, 100}, point{110, 110}},
    };
    auto g = build_grid(boxes);

    // Query first
    int count = 0;
    g.query_intersects(box{point{5, 5}, point{15, 15}}, [&](size_t) { count++; });
    EXPECT_EQ(count, 1);

    // Query middle (empty)
    count = 0;
    g.query_intersects(box{point{50, 50}, point{60, 60}}, [&](size_t) { count++; });
    EXPECT_EQ(count, 0);
}

TEST(UniformGrid, TwoOverlapping) {
    std::vector<box> boxes = {
        {point{0, 0}, point{10, 10}},
        {point{5, 5}, point{15, 15}},
    };
    auto g = build_grid(boxes);

    // Query overlapping region
    int count = 0;
    g.query_intersects(box{point{8, 8}, point{12, 12}}, [&](size_t) { count++; });
    EXPECT_EQ(count, 2);
}

TEST(UniformGrid, LargeGridThousands) {
    // 2000 items — should trigger grid mode
    std::vector<std::pair<box, size_t>> items;
    items.reserve(2000);
    for (size_t i = 0; i < 2000; i++) {
        int32_t x = (int32_t)(2 * i);
        items.push_back({{point{x, x}, point{x + 2, x + 2}}, i});
    }
    grid<size_t> g(std::move(items), 0);

    // Query overlapping a middle pair
    int count = 0;
    g.query_intersects(box{point{5, 5}, point{6, 6}}, [&](size_t) { count++; });
    EXPECT_GE(count, 1);
}

TEST(UniformGrid, LargeGridWithAdjacentPairs) {
    // Overlapping pairs: (0,0)-(2,2) and (1,1)-(3,3), (2,2)-(4,4) and (3,3)-(5,5), ...
    std::vector<std::pair<box, size_t>> items;
    size_t n = 1000;
    items.reserve(n * 2);
    for (size_t i = 0; i < n; i++) {
        items.push_back({{point{0 + 2 * (int32_t)i, 0 + 2 * (int32_t)i},
                          point{2 + 2 * (int32_t)i, 2 + 2 * (int32_t)i}}, i * 2});
        items.push_back({{point{1 + 2 * (int32_t)i, 1 + 2 * (int32_t)i},
                          point{3 + 2 * (int32_t)i, 3 + 2 * (int32_t)i}}, i * 2 + 1});
    }
    grid<size_t> g(std::move(items), 0);

    // Query each overlapping pair
    for (size_t i = 0; i < n; i++) {
        int count = 0;
        g.query_intersects(
            box{point{1 + 2 * (int32_t)i, 1 + 2 * (int32_t)i},
                point{2 + 2 * (int32_t)i, 2 + 2 * (int32_t)i}},
            [&](size_t) { count++; });
        EXPECT_GE(count, 2);  // at least both squares in the pair
    }
}

TEST(UniformGrid, PointQueryTouchingCorner) {
    // Box at (10,10)-(20,20), query at exact corner (20,20)
    std::vector<box> boxes = {{point{10, 10}, point{20, 20}}};
    auto g = build_grid(boxes);

    int count = 0;
    g.query_intersects(box{point{20, 20}, point{20, 20}}, [&](size_t) { count++; });
    EXPECT_EQ(count, 1);
}

TEST(UniformGrid, PointQueryOutside) {
    std::vector<box> boxes = {{point{10, 10}, point{20, 20}}};
    auto g = build_grid(boxes);

    int count = 0;
    g.query_intersects(box{point{21, 21}, point{21, 21}}, [&](size_t) { count++; });
    EXPECT_EQ(count, 0);
}

TEST(UniformGrid, QueryLargerThanData) {
    // Query box that covers entire data range
    std::vector<box> boxes = {
        {point{0, 0}, point{10, 10}},
        {point{20, 20}, point{30, 30}},
        {point{40, 40}, point{50, 50}},
    };
    auto g = build_grid(boxes);

    int count = 0;
    g.query_intersects(box{point{0, 0}, point{100, 100}}, [&](size_t) { count++; });
    EXPECT_GE(count, 3);
}

TEST(UniformGrid, HighCoordinateValues) {
    // Coordinates near INT32_MAX (after biasing, these are near the uint32 mid-range)
    int32_t base = 0;
    std::vector<box> boxes = {
        {point{base + 0, base + 0}, point{base + 10, base + 10}},
        {point{base + 5, base + 5}, point{base + 15, base + 15}},
        {point{base + 20, base + 20}, point{base + 30, base + 30}},
    };
    auto g = build_grid(boxes);

    // Query overlapping pair
    int count = 0;
    g.query_intersects(box{point{base + 7, base + 7}, point{base + 12, base + 12}},
                       [&](size_t) { count++; });
    EXPECT_EQ(count, 2);
}

TEST(UniformGrid, QuerySmallBoxReturnsSelf) {
    // A single tiny box, query with the same box — should return it
    std::vector<box> boxes = {{point{5, 5}, point{6, 6}}};
    auto g = build_grid(boxes);

    int count = 0;
    g.query_intersects(boxes[0], [&](size_t) { count++; });
    EXPECT_GE(count, 1);
}

TEST(UniformGrid, DegenerateBoxSinglePoint) {
    // Box with zero area (point at (42, 42))
    std::vector<box> boxes = {{point{42, 42}, point{42, 42}}};
    auto g = build_grid(boxes);

    // Query exactly that point
    int count = 0;
    g.query_intersects(box{point{42, 42}, point{42, 42}}, [&](size_t) { count++; });
    EXPECT_EQ(count, 1);
}
