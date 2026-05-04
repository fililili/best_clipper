#include <gtest/gtest.h>
#include "uint32_adaptor.hpp"

// ============================================================================
// Helpers — verify against known WKT
// ============================================================================

void test_union(const std::string& a_wkt, const std::string& b_wkt,
                const std::string& expected_wkt) {
    multi_polygon_s32 a, b, expected;
    bg::read_wkt(a_wkt, a);     ASSERT_TRUE(bg::is_valid(a));
    bg::read_wkt(b_wkt, b);     ASSERT_TRUE(bg::is_valid(b));
    bg::read_wkt(expected_wkt, expected);
    ASSERT_TRUE(bg::is_valid(expected));
    auto result = add(a, b);
    EXPECT_TRUE(bg::is_valid(result));
    EXPECT_TRUE(bg::equals(result, expected))
        << "Result:   " << bg::wkt(result) << "\n"
        << "Expected: " << bg::wkt(expected);
}

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

// ============================================================================
// Union tests
// ============================================================================

TEST(BasicTest, Union) {
    test_union(
        "MULTIPOLYGON(((-59 867,-36 492,-182 486,-59 867)))",
        "MULTIPOLYGON(((-220 877,-54 821,-402 541,-808 638,-220 877)))",
        "MULTIPOLYGON(((-220 877,-72 827,-59 867,-56 822,-54 821,-56 819,-36 492,-182 486,-81 799,-402 541,-808 638,-220 877)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))",
        "MULTIPOLYGON(((2 2, 2 4, 4 4, 4 2, 2 2)))",
        "MULTIPOLYGON(((0 0, 0 3, 2 3, 2 4, 4 4, 4 2, 3 2, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 2, 5 1, 5 0, 0 0)))",
        "MULTIPOLYGON(((3 1, 0 2, 5 1, 3 1)))",
        "MULTIPOLYGON(((0 2,3 1,5 1,5 0,0 0,0 2)))"
    );
    test_union(
        "MULTIPOLYGON(((-1 -1, -1 3, 3 3, 3 -1, -1 -1), (0 0, 2 0, 2 2, 0 2, 0 0)))",
        "MULTIPOLYGON(((1 1, 1 4, 4 4, 4 1, 1 1)))",
        "MULTIPOLYGON(((-1 3,1 3,1 4,4 4,4 1,3 1,3 -1,-1 -1,-1 3),(2 0,2 1,1 1,1 2,0 2,0 0,2 0)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1), (6 6, 8 6, 8 8, 6 8, 6 6)))",
        "MULTIPOLYGON(((2 2, 2 7, 7 7, 7 2, 2 2)))",
        "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 2, 2 2, 2 3, 1 3, 1 1), (8 8, 6 8, 6 7, 7 7, 7 6, 8 6, 8 8)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 1 1, 2 1, 2 2, 3 3, 3 0, 0 0)))",
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 2 2, 1 2, 1 1, 0 0)))",
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))"
    );
    test_union(
        "MULTIPOLYGON(((-1461 -786,-1417 -833,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1030 -1246,-1608 -939,-1461 -786),(-1058 -1133,-1244 -963,-1301 -1039,-1058 -1133)),((-578 -670,-688 -679,-802 -443,-826 -400,-578 -670)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))",
        "MULTIPOLYGON(((-1450 -1280, -1450 -800, -1200 -1000, -1000 -1280, -1450 -1280)))",
        "MULTIPOLYGON(((-1461 -786,-1442 -807,-1410 -832,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1026 -1244,-1000 -1280,-1450 -1280,-1450 -1023,-1608 -939,-1461 -786),(-1245 -964,-1228 -977,-1244 -963,-1245 -964),(-1123 -1108,-1058 -1133,-1193 -1009,-1123 -1108)),((-826 -400,-578 -670,-688 -679,-802 -443,-826 -400)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))"
    );
}

// ============================================================================
// Rectangle union / self_or — many overlapping squares
// ============================================================================

void test_union_rectangle(int size) {
    multi_polygon_s32 a, b;
    for (int i = 0; i < size; i++) {
        a.emplace_back(polygon_s32{ {{0 + 2 * i, 0 + 2 * i},
                                 {0 + 2 * i, 2 + 2 * i},
                                 {2 + 2 * i, 2 + 2 * i},
                                 {2 + 2 * i, 0 + 2 * i},
                                 {0 + 2 * i, 0 + 2 * i}} });
        b.emplace_back(polygon_s32{ {{1 + 2 * i, 1 + 2 * i},
                                 {1 + 2 * i, 3 + 2 * i},
                                 {3 + 2 * i, 3 + 2 * i},
                                 {3 + 2 * i, 1 + 2 * i},
                                 {1 + 2 * i, 1 + 2 * i}} });
    }
    EXPECT_TRUE(bg::is_valid(a));
    EXPECT_TRUE(bg::is_valid(b));
    EXPECT_TRUE(bg::equals(self_or(a), a));
    EXPECT_TRUE(bg::equals(self_or(b), b));
    auto result = add(a, b);
    EXPECT_TRUE(bg::is_valid(result));
    EXPECT_DOUBLE_EQ(bg::area(result), 1 + 6 * size);
}

void test_self_or_rectangle(int size) {
    multi_polygon_s32 poly;
    for (int i = 0; i < size; i++) {
        poly.emplace_back(polygon_s32{ {{0 + i, 0 + i},
                                    {0 + i, 2 + i},
                                    {2 + i, 2 + i},
                                    {2 + i, 0 + i},
                                    {0 + i, 0 + i}} });
    }
    auto result = self_or(poly);
    EXPECT_TRUE(bg::is_valid(result));
    EXPECT_DOUBLE_EQ(bg::area(result), 1 + 3 * size);
}

TEST(BasicUnion, RectangleUnion) {
    test_union_rectangle(100);
    test_union_rectangle(3);
    test_union_rectangle(1010);
    test_union_rectangle(300);
    test_union_rectangle(2521);
}

TEST(BasicUnion, RectangleSelfOr) {
    test_self_or_rectangle(100);
    test_self_or_rectangle(3);
    test_self_or_rectangle(1010);
    test_self_or_rectangle(300);
    test_self_or_rectangle(2521);
}

// ============================================================================
// Intersection tests
// ============================================================================

TEST(BasicTest, Intersection) {
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

// ============================================================================
// Difference tests
// ============================================================================

TEST(BasicTest, Difference) {
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

// ============================================================================
// connected_components unit tests
// ============================================================================

TEST(ConnectedComponents, EmptyGraph) {
    auto comp = connected_components(5, {});
    EXPECT_EQ(comp.size(), 5u);
    for (std::size_t i = 0; i < 5; i++)
        EXPECT_EQ(comp[i], i) << "vertex " << i;
}

TEST(ConnectedComponents, SingleEdge) {
    auto comp = connected_components(4, {{0, 1}});
    EXPECT_EQ(comp[0], comp[1]);
    EXPECT_NE(comp[0], comp[2]);
    EXPECT_NE(comp[0], comp[3]);
    EXPECT_NE(comp[2], comp[3]);
}

TEST(ConnectedComponents, Chain) {
    auto comp = connected_components(5, {{0, 1}, {1, 2}, {2, 3}});
    EXPECT_EQ(comp[0], comp[1]);
    EXPECT_EQ(comp[0], comp[2]);
    EXPECT_EQ(comp[0], comp[3]);
    EXPECT_NE(comp[0], comp[4]);
}

TEST(ConnectedComponents, Disconnected) {
    auto comp = connected_components(6, {{0, 1}, {2, 3}, {4, 5}});
    EXPECT_EQ(comp[0], comp[1]);
    EXPECT_EQ(comp[2], comp[3]);
    EXPECT_EQ(comp[4], comp[5]);
    EXPECT_NE(comp[0], comp[2]);
    EXPECT_NE(comp[0], comp[4]);
    EXPECT_NE(comp[2], comp[4]);
}

TEST(ConnectedComponents, DuplicateEdges) {
    auto comp = connected_components(3, {{0, 1}, {0, 1}, {1, 0}});
    EXPECT_EQ(comp[0], comp[1]);
    EXPECT_NE(comp[0], comp[2]);
}

TEST(ConnectedComponents, SelfLoop) {
    auto comp = connected_components(3, {{0, 0}, {0, 1}});
    EXPECT_EQ(comp[0], comp[1]);
    EXPECT_NE(comp[0], comp[2]);
}

TEST(ConnectedComponents, StarGraph) {
    auto comp = connected_components(5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}});
    for (std::size_t i = 1; i < 5; i++)
        EXPECT_EQ(comp[0], comp[i]);
}

TEST(ConnectedComponents, IsolatedVertex) {
    auto comp = connected_components(4, {{0, 1}});
    EXPECT_EQ(comp[0], comp[1]);
    EXPECT_NE(comp[2], comp[0]);
    EXPECT_NE(comp[3], comp[0]);
    EXPECT_NE(comp[2], comp[3]);
}

TEST(ConnectedComponents, LargeChain) {
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    for (std::size_t i = 0; i + 1 < 1000; i++)
        edges.emplace_back(i, i + 1);
    auto comp = connected_components(1000, edges);
    for (std::size_t i = 1; i < 1000; i++)
        EXPECT_EQ(comp[0], comp[i]);
    EXPECT_EQ(comp[0], 0u);
}
