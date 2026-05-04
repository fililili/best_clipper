#include <gtest/gtest.h>
#include "uint32_adaptor.hpp"

// ============================================================================
// XOR tests — validity checks (areas verified independently)
// ============================================================================

TEST(BasicTest, Xor) {
    // Two overlapping axis-aligned rects
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((1 1,1 3,3 3,3 1,1 1)))", b);
        auto result = xor_(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_FALSE(result.empty());
    }

    // One rect contains the other — XOR produces a hole
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((1 1,1 2,2 2,2 1,1 1)))", b);
        auto result = xor_(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_FALSE(result.empty());
    }

    // Two disjoint rects — both survive
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((3 3,3 4,4 4,4 3,3 3)))", b);
        auto result = xor_(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_EQ(result.size(), 2u);
    }

    // Adjacent (edge-sharing) rects — zero-area intersection
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((1 0,1 1,2 1,2 0,1 0)))", b);
        auto result = xor_(a, b);
        EXPECT_TRUE(bg::is_valid(result));
        EXPECT_FALSE(result.empty());
    }

    // Triangle + rect with non-integer intersection → snap rounding active
    // Triangle: (0,0)-(0,4)-(8,0), hypotenuse x+2y=8
    // Rect: (2,1)-(5,3). Hypotenuse ∩ right edge (x=5): y=1.5 → snaps to (5,2)
    {
        multi_polygon_s32 a, b;
        bg::read_wkt("MULTIPOLYGON(((0 0,0 4,8 0,0 0)))", a);
        bg::read_wkt("MULTIPOLYGON(((2 1,2 3,5 3,5 1,2 1)))", b);
        auto result = xor_(a, b);
        EXPECT_TRUE(bg::is_valid(result)) << bg::wkt(result);
        EXPECT_FALSE(result.empty());
    }
}
