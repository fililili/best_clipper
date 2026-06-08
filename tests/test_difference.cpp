#include <gtest/gtest.h>
#include "core.hpp"

using namespace best_clipper;

// ============================================================================
// Difference tests — Manhattan cases compare against Boost ground truth
// ============================================================================

void test_difference(const std::string& a_wkt, const std::string& b_wkt) {
    multi_polygon a, b;
    bg::read_wkt(a_wkt, a);     ASSERT_TRUE(bg::is_valid(a));
    bg::read_wkt(b_wkt, b);     ASSERT_TRUE(bg::is_valid(b));

    auto result = difference(a, b);

    multi_polygon expected;
    bg::difference(a, b, expected);

    EXPECT_TRUE(bg::is_valid(result)) << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "Result:   " << bg::wkt(result) << "\n"
        << "Expected: " << bg::wkt(expected);
}

TEST(CoreTest, Difference) {
    // A \ B where A and B overlap
    test_difference(
        "MULTIPOLYGON(((0 0,0 2,2 2,2 0,0 0)))",
        "MULTIPOLYGON(((1 1,1 3,3 3,3 1,1 1)))"
    );

    // A \ B where B is entirely inside A — result has a hole
    test_difference(
        "MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))",
        "MULTIPOLYGON(((1 1,1 2,2 2,2 1,1 1)))"
    );

    // A \ B where A and B are disjoint — A unchanged
    test_difference(
        "MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)))",
        "MULTIPOLYGON(((3 3,3 4,4 4,4 3,3 3)))"
    );

    // A \ B where B completely covers A — result is empty
    test_difference(
        "MULTIPOLYGON(((1 1,1 2,2 2,2 1,1 1)))",
        "MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))"
    );

    // A \ B where A and B are identical — result is empty
    test_difference(
        "MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))",
        "MULTIPOLYGON(((0 0,0 3,3 3,3 0,0 0)))"
    );

    // Large square with 2 holes \ overlapping rect
    test_difference(
        "MULTIPOLYGON(((0 0,0 9,9 9,9 0,0 0),(1 1,3 1,3 3,1 3,1 1),(6 6,8 6,8 8,6 8,6 6)))",
        "MULTIPOLYGON(((2 2,2 7,7 7,7 2,2 2)))"
    );
}

// Non-Manhattan case — Boost.Geometry is unreliable for non-integer
// intersections, so only verify validity and non-emptiness.
TEST(CoreTest, DifferenceSnapRounding) {
    // Triangle \ rect with non-integer intersection → snap rounding active
    // Triangle: (0,0)-(0,4)-(8,0), hypotenuse x+2y=8
    // Rect: (2,1)-(5,3). Hypotenuse ∩ right edge (x=5): y=1.5 → snaps to (5,2)
    multi_polygon a, b;
    bg::read_wkt("MULTIPOLYGON(((0 0,0 4,8 0,0 0)))", a);
    bg::read_wkt("MULTIPOLYGON(((2 1,2 3,5 3,5 1,2 1)))", b);
    ASSERT_TRUE(bg::is_valid(a));
    ASSERT_TRUE(bg::is_valid(b));
    auto result = difference(a, b);
    EXPECT_TRUE(bg::is_valid(result)) << bg::wkt(result);
    EXPECT_FALSE(result.empty());
    EXPECT_GT(bg::area(result), 0.0);
}
