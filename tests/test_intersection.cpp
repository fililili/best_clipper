#include <gtest/gtest.h>
#include "uint32_adaptor.hpp"

// ============================================================================
// Intersection tests — verify against known WKT
// ============================================================================

void test_intersection(const std::string& a_wkt, const std::string& b_wkt,
                       const std::string& expected_wkt) {
    multi_polygon_s32 a, b, expected;
    bg::read_wkt(a_wkt, a);     ASSERT_TRUE(bg::is_valid(a));
    bg::read_wkt(b_wkt, b);     ASSERT_TRUE(bg::is_valid(b));
    bg::read_wkt(expected_wkt, expected);
    ASSERT_TRUE(bg::is_valid(expected));
    auto result = intersection(a, b);
    EXPECT_TRUE(bg::is_valid(result));
    EXPECT_TRUE(bg::equals(result, expected))
        << "Result:   " << bg::wkt(result) << "\n"
        << "Expected: " << bg::wkt(expected);
}

TEST(CoreTest, Intersection) {
    test_intersection(
        "MULTIPOLYGON(((0 0, 0 2, 2 2, 2 0, 0 0)))",
        "MULTIPOLYGON(((1 1, 1 3, 3 3, 3 1, 1 1)))",
        "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))"
    );
    test_intersection(
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0)))",
        "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))",
        "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))"
    );
    test_intersection(
        "MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))",
        "MULTIPOLYGON(((2 2, 2 3, 3 3, 3 2, 2 2)))",
        "MULTIPOLYGON()"
    );
    test_intersection(
        "MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))",
        "MULTIPOLYGON(((1 0, 1 1, 2 1, 2 0, 1 0)))",
        "MULTIPOLYGON()"
    );
    test_intersection(
        "MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))",
        "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))",
        "MULTIPOLYGON()"
    );
    test_intersection(
        "MULTIPOLYGON(((0 0, 0 10, 10 10, 10 0, 0 0)), ((20 0, 20 10, 30 10, 30 0, 20 0)))",
        "MULTIPOLYGON(((5 5, 5 15, 15 15, 15 5, 5 5)))",
        "MULTIPOLYGON(((5 5, 5 10, 10 10, 10 5, 5 5)))"
    );
}
