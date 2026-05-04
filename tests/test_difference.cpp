#include <gtest/gtest.h>
#include "uint32_adaptor.hpp"

// ============================================================================
// Difference tests
// ============================================================================

TEST(CoreTest, Difference) {
    // A \ B where A and B overlap
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((1 1,1 3,3 3,3 1,1 1)))", b);
        auto result = difference(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_FALSE(result.empty());
    }

    // A \ B where B is entirely inside A — result has a hole
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((1 1,1 2,2 2,2 1,1 1)))", b);
        auto result = difference(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_FALSE(result.empty());
    }

    // A \ B where A and B are disjoint — A unchanged
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((3 3,3 4,4 4,4 3,3 3)))", b);
        auto result = difference(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_FALSE(result.empty());
    }

    // A \ B where B completely covers A — result is empty
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((1 1,1 2,2 2,2 1,1 1)))", a);
        bg::read_wkt("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))", b);
        auto result = difference(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_TRUE(result.empty());
    }

    // A \ B where A and B are identical — result is empty
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))", b);
        auto result = difference(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_TRUE(result.empty());
    }

    // Triangle \ rect with non-integer intersection → snap rounding active
    // Triangle: (0,0)-(0,4)-(8,0), hypotenuse x+2y=8
    // Rect: (2,1)-(5,3). Hypotenuse ∩ right edge (x=5): y=1.5 → snaps to (5,2)
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 4,8 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((2 1,2 3,5 3,5 1,2 1)))", b);
        auto result = difference(a, b);
        EXPECT_TRUE(bg::is_valid(result)) << bg::wkt(result);
        EXPECT_FALSE(result.empty());
    }
}
