#include <gtest/gtest.h>
#include "uint32_adaptor.hpp"

// All polygons below are OGC valid: no self-intersections, no overlapping
// components, correct ring orientation. All edges are axis-aligned with integer
// coordinates, so every intersection already falls on an integer grid point —
// snap rounding is a no-op for these inputs.
//
// Because the input is OGC valid, Boost.Geometry is reliable ground truth.
// Boost is only unreliable when the input is invalid (e.g. self-overlapping
// multi_polygon components). For valid input, Boost produces correct results.

// Use Boost.Geometry's union_ and intersection as ground truth.
multi_polygon_s32 bg_self_union(const multi_polygon_s32& p) {
    if (p.empty()) return p;
    multi_polygon_s32 out;
    out.push_back(p[0]);
    for (size_t i = 1; i < p.size(); i++) {
        multi_polygon_s32 tmp;
        bg::union_(out, multi_polygon_s32{p[i]}, tmp);
        out = std::move(tmp);
    }
    bg::correct(out);
    return out;
}

multi_polygon_s32 bg_union(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    auto au = bg_self_union(a);
    auto bu = bg_self_union(b);
    multi_polygon_s32 out;
    bg::union_(au, bu, out);
    bg::correct(out);
    return out;
}

multi_polygon_s32 bg_intersection(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    auto au = bg_self_union(a);
    auto bu = bg_self_union(b);
    multi_polygon_s32 out;
    bg::intersection(au, bu, out);
    bg::correct(out);
    return out;
}

multi_polygon_s32 bg_xor(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    auto au = bg_self_union(a);
    auto bu = bg_self_union(b);
    multi_polygon_s32 un, inter, out;
    bg::union_(au, bu, un);
    bg::intersection(au, bu, inter);
    bg::difference(un, inter, out);
    bg::correct(out);
    return out;
}

multi_polygon_s32 bg_difference(const multi_polygon_s32& a, const multi_polygon_s32& b) {
    auto au = bg_self_union(a);
    auto bu = bg_self_union(b);
    multi_polygon_s32 out;
    bg::difference(au, bu, out);
    bg::correct(out);
    return out;
}

// Helper: create a rectangle ring. CCW order.
ring_s32 make_rect(int32_t x1, int32_t y1, int32_t x2, int32_t y2) {
    ring_s32 r;
    r.push_back(point_s32{x1, y1});
    r.push_back(point_s32{x1, y2});
    r.push_back(point_s32{x2, y2});
    r.push_back(point_s32{x2, y1});
    r.push_back(point_s32{x1, y1}); // close
    return r;
}

// Helper: create a CW ring (hole)
ring_s32 make_hole(int32_t x1, int32_t y1, int32_t x2, int32_t y2) {
    ring_s32 r;
    r.push_back(point_s32{x1, y1});
    r.push_back(point_s32{x2, y1});
    r.push_back(point_s32{x2, y2});
    r.push_back(point_s32{x1, y2});
    r.push_back(point_s32{x1, y1}); // close
    return r;
}

polygon_s32 make_rect_poly(int32_t x1, int32_t y1, int32_t x2, int32_t y2) {
    polygon_s32 p;
    p.outer() = make_rect(x1, y1, x2, y2);
    return p;
}

// ============================================================================
// Manhattan Union Tests
// ============================================================================

TEST(ManhattanUnion, TwoSeparateRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(5, 5, 7, 7));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, OverlappingRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 3, 3));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, AdjacentRects) {
    // Sharing an edge
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 0, 4, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 8.0);
}

TEST(ManhattanUnion, CornerTouching) {
    // Sharing a corner only
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 4, 4));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, OneContainsOther) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 10, 10));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 100.0);
}

TEST(ManhattanUnion, RectWithHole) {
    // Outer rect has a hole; union with a rect that overlaps the hole
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 6, 6);
    p1.inners().push_back(make_hole(2, 2, 4, 4));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(3, 3, 8, 8));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, RectPartiallyFillsHole) {
    // Outer rect with hole; union with rect that partially covers the hole
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 8, 8);
    p1.inners().push_back(make_hole(2, 2, 6, 6));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(4, 4, 7, 7));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, TwoHolesOneFilled) {
    // Outer rect with 2 holes; union with rect that fills one hole
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 10, 10);
    p1.inners().push_back(make_hole(1, 1, 4, 4));
    p1.inners().push_back(make_hole(6, 1, 9, 4));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 8, 3));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, NestedMultiPolygon) {
    // Two separate polygons in first, one overlapping rect in second
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 3, 3));
    a.push_back(make_rect_poly(5, 0, 8, 3));
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 7, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 20.0); // 9+9+2(gap fill) = 20
}

TEST(ManhattanUnion, GridPattern) {
    // 3x3 grid of squares with a diagonal cross
    multi_polygon_s32 a;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a.push_back(make_rect_poly(i * 10, j * 10, i * 10 + 4, j * 10 + 4));
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(2, 2, 22, 22));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, LShaped) {
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 5, 2);
        polygon_s32 q; q.outer() = make_rect(0, 2, 2, 5);
        a.push_back(p); a.push_back(q);
        a = bg_self_union(a); // make it a valid L-shape
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(3, 0, 6, 3));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, CrossShape) {
    // A cross shape (like a plus sign)
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(3, 0, 5, 10));
    a.push_back(make_rect_poly(0, 3, 8, 5));
    a = bg_self_union(a);
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(2, 2, 6, 6));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanUnion, ManyOverlappingSquares) {
    // 5 overlapping squares forming a staircase
    multi_polygon_s32 a;
    for (int i = 0; i < 5; i++)
        a.push_back(make_rect_poly(i * 2, i * 2, i * 2 + 3, i * 2 + 3));
    multi_polygon_s32 b;
    for (int i = 0; i < 5; i++)
        b.push_back(make_rect_poly(i * 2 + 1, i * 2 - 1, i * 2 + 4, i * 2 + 2));
    // a and b each contain internally-overlapping squares — multi_polygon
    // validity requires non-overlapping polygons, skip input validity check
    ASSERT_TRUE(bg::is_valid(bg_self_union(a))) << "self-union a invalid";
    ASSERT_TRUE(bg::is_valid(bg_self_union(b))) << "self-union b invalid";

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

// ============================================================================
// Manhattan Intersection Tests
// ============================================================================

TEST(ManhattanIntersection, TwoSeparateRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(5, 5, 7, 7));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_TRUE(result.empty());
}

TEST(ManhattanIntersection, OverlappingRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 3, 3));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 1.0);
}

TEST(ManhattanIntersection, AdjacentRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 0, 4, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanIntersection, CornerTouching) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 4, 4));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanIntersection, OneContainsOther) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 10, 10));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 9.0);
}

TEST(ManhattanIntersection, RectWithHole) {
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 6, 6);
    p1.inners().push_back(make_hole(2, 2, 4, 4));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(1, 1, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanIntersection, MultiPolygonBothSides) {
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 3, 3));
    a.push_back(make_rect_poly(5, 0, 8, 3));
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 7, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 4.0); // (2x1)+(2x1) = 4
}

// ============================================================================
// Self-union Tests (Manhattan)
// ============================================================================

TEST(ManhattanSelfOr, TwoOverlappingRects) {
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 3, 3));
    a.push_back(make_rect_poly(2, 2, 5, 5));
    // a has overlapping parts – input itself may not be valid as multi_polygon
    // because bg considers overlapping polygons in multi_polygon invalid

    auto result = self_or(a);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_DOUBLE_EQ(bg::area(result), 17.0); // 9+9-1(overlap) = 17
}

TEST(ManhattanSelfOr, DisjointRects) {
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 2, 2));
    a.push_back(make_rect_poly(4, 4, 6, 6));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);

    auto result = self_or(a);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_DOUBLE_EQ(bg::area(result), 8.0);
    EXPECT_EQ(result.size(), 2u);
}

TEST(ManhattanSelfOr, NestedRect) {
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 10, 10));
    a.push_back(make_rect_poly(2, 2, 5, 5));
    // overlapping — input invalid as multi_polygon per bg convention

    auto result = self_or(a);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_DOUBLE_EQ(bg::area(result), 100.0);
    EXPECT_EQ(result.size(), 1u);
}

TEST(ManhattanSelfOr, ThreeOverlapping) {
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 3, 3));
    a.push_back(make_rect_poly(2, 1, 5, 4));
    a.push_back(make_rect_poly(1, 2, 4, 6));
    // overlapping — input invalid as multi_polygon per bg convention

    auto result = self_or(a);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
}

TEST(ManhattanSelfOr, LShapedSelfOr) {
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 5, 2));
    a.push_back(make_rect_poly(0, 2, 2, 5));
    // edge-touching inputs — bg considers overlapping/touching multi_polygon invalid
    ASSERT_TRUE(bg::is_valid(bg_self_union(a))) << "self-union a invalid";

    auto result = self_or(a);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_EQ(result.size(), 1u);
    EXPECT_DOUBLE_EQ(bg::area(result), 16.0); // 10+6 (line overlap only, no area loss)
}

// ============================================================================
// Edge cases
// ============================================================================

TEST(ManhattanEdgeCase, ZeroAreaIntersection) {
    // Intersection is a line segment (zero area)
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 0, 4, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    // Should be empty (zero-area intersection along shared edge)
    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanEdgeCase, VeryNarrowOverlap) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 10, 10));
    multi_polygon_s32 b; b.push_back(make_rect_poly(9, 0, 11, 10));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanEdgeCase, ExactlyFittingHole) {
    // Union fills the hole exactly
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 10, 10);
    p1.inners().push_back(make_hole(2, 2, 6, 6));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 6, 6));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 100.0);
    EXPECT_EQ(result.size(), 1u);
}

TEST(ManhattanEdgeCase, LastColumnOverlap) {
    // Two rects that overlap only in the last column of pixels
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 5, 5));
    multi_polygon_s32 b; b.push_back(make_rect_poly(4, 2, 8, 3));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanEdgeCase, MultipleHoles) {
    polygon_s32 p; p.outer() = make_rect(0, 0, 10, 10);
    p.inners().push_back(make_hole(1, 1, 3, 3));
    p.inners().push_back(make_hole(4, 1, 6, 3));
    p.inners().push_back(make_hole(7, 1, 9, 3));
    p.inners().push_back(make_hole(1, 5, 3, 9));
    p.inners().push_back(make_hole(5, 6, 9, 7));
    multi_polygon_s32 a; a.push_back(p);
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(2, 2, 8, 8));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanEdgeCase, SharingMultipleEdges) {
    // Two L-shapes sharing a partial boundary
    multi_polygon_s32 a;
    a.push_back(make_rect_poly(0, 0, 4, 2));
    a.push_back(make_rect_poly(0, 2, 2, 4));
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(2, 0, 6, 2));
    b.push_back(make_rect_poly(4, 2, 6, 4));
    // edge-touching inputs — bg considers touching multi_polygon invalid
    ASSERT_TRUE(bg::is_valid(bg_self_union(a))) << "self-union a invalid";
    ASSERT_TRUE(bg::is_valid(bg_self_union(b))) << "self-union b invalid";

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

// ============================================================================
// Large Scale Tests
// ============================================================================

TEST(ManhattanLarge, RectWith100Holes) {
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 10; i++)
            for (int j = 0; j < 10; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 10, 8 + i * 10, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWith20Holes) {
    // Scaled-down version of RectWith100Holes for debugging
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 5; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 20, 8 + i * 10, 8 + j * 20));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithCoincidentHole) {
    // Hole left edge at x=52 coincides with B's right edge — the root issue.
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        p.inners().push_back(make_hole(52, 2, 58, 8)); // hole column 5, left at x=52
        p.inners().push_back(make_hole(2, 2, 8, 8));   // hole column 0, inside B
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99)); // B's right edge at x=52 coincides with hole left
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithTwoCoincidentHoles) {
    // Two holes sharing left boundary with B at x=52
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        p.inners().push_back(make_hole(52, 2, 58, 8));   // hole 0, x=52..58, y=2..8
        p.inners().push_back(make_hole(52, 12, 58, 18)); // hole 1, x=52..58, y=12..18
        p.inners().push_back(make_hole(2, 2, 8, 8));     // inside B
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithTenCoincidentHoles) {
    // 10 holes all with left edge at x=52 (coincident with B)
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(52, 2 + j * 10, 58, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithThreeColumns) {
    // Column 4 (inside B), column 5 (coincident), column 6 (outside B)
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        // Column 4: inside B, should be filled
        for (int j = 0; j < 3; j++)
            p.inners().push_back(make_hole(42, 2 + j * 10, 48, 8 + j * 10));
        // Column 5: coincident with B, should NOT be filled
        for (int j = 0; j < 3; j++)
            p.inners().push_back(make_hole(52, 2 + j * 10, 58, 8 + j * 10));
        // Column 6: outside B, should NOT be filled
        for (int j = 0; j < 3; j++)
            p.inners().push_back(make_hole(62, 2 + j * 10, 68, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithSixColumnsPlusColumn9) {
    // Columns 0-5 + column 9: does adding a far-away column trigger the bug?
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 10; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 10, 8 + i * 10, 8 + j * 10));
        // Add column 9 (far right)
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(92, 2 + j * 10, 98, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithSixColumns) {
    // Columns 0-5, 10 rows each = 60 holes. Cols 0-4 inside B, col 5 coincident.
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 10; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 10, 8 + i * 10, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithSixColumnsPlusColumn9Intersection) {
    // Intersection: columns 0-5 + column 9 ∩ rect at x=1..52
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 10; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 10, 8 + i * 10, 8 + j * 10));
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(92, 2 + j * 10, 98, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = intersection(a, b);
    auto expected = bg_intersection(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithTwoColumnsFarApart) {
    // Only columns 0 and 9 — max distance between hole groups
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(2, 2 + j * 10, 8, 8 + j * 10));
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(92, 2 + j * 10, 98, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 10, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

TEST(ManhattanLarge, RectWithSixColumnsPlusColumn8And9) {
    // Two far-away columns (8 and 9) — multiple isolated groups
    multi_polygon_s32 a;
    {
        polygon_s32 p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 10; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 10, 8 + i * 10, 8 + j * 10));
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(82, 2 + j * 10, 88, 8 + j * 10));
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(92, 2 + j * 10, 98, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon_s32 b;
    b.push_back(make_rect_poly(1, 1, 52, 99));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = add(a, b);
    auto expected = bg_union(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::is_valid(expected));
    EXPECT_DOUBLE_EQ(bg::area(result), bg::area(expected));
}

// ============================================================================
// Manhattan XOR Tests
// ============================================================================

TEST(ManhattanXor, TwoSeparateRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(5, 5, 7, 7));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = xor_(a, b);
    auto expected = bg_xor(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 8.0);
}

TEST(ManhattanXor, OverlappingRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 3, 3));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = xor_(a, b);
    auto expected = bg_xor(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 16.0);
}

TEST(ManhattanXor, AdjacentRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 0, 4, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = xor_(a, b);
    auto expected = bg_xor(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanXor, OneContainsOther) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 10, 10));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = xor_(a, b);
    auto expected = bg_xor(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 91.0);
}

TEST(ManhattanXor, RectWithHole) {
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 6, 6);
    p1.inners().push_back(make_hole(2, 2, 4, 4));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(3, 3, 8, 8));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = xor_(a, b);
    auto expected = bg_xor(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

// ============================================================================
// Manhattan Difference Tests
// ============================================================================

TEST(ManhattanDifference, TwoSeparateRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(5, 5, 7, 7));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 4.0);
}

TEST(ManhattanDifference, OverlappingRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 3, 3));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 8.0);
}

TEST(ManhattanDifference, AdjacentRects) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 2, 2));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 0, 4, 2));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 4.0);
}

TEST(ManhattanDifference, OneContainsOther) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(0, 0, 10, 10));
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 5, 5));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 91.0);
}

TEST(ManhattanDifference, BcontainsA) {
    multi_polygon_s32 a; a.push_back(make_rect_poly(2, 2, 5, 5));
    multi_polygon_s32 b; b.push_back(make_rect_poly(0, 0, 10, 10));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_TRUE(result.empty());
}

TEST(ManhattanDifference, RectWithHole) {
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 6, 6);
    p1.inners().push_back(make_hole(2, 2, 4, 4));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(3, 3, 8, 8));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
}

TEST(ManhattanDifference, RectMinusHoleFiller) {
    polygon_s32 p1; p1.outer() = make_rect(0, 0, 10, 10);
    p1.inners().push_back(make_hole(2, 2, 6, 6));
    multi_polygon_s32 a; a.push_back(p1);
    multi_polygon_s32 b; b.push_back(make_rect_poly(2, 2, 6, 6));
    ASSERT_TRUE(bg::is_valid(a)) << "a: " << bg::wkt(a);
    ASSERT_TRUE(bg::is_valid(b)) << "b: " << bg::wkt(b);

    auto result = difference(a, b);
    auto expected = bg_difference(a, b);

    EXPECT_TRUE(bg::is_valid(result)) << "result: " << bg::wkt(result);
    EXPECT_TRUE(bg::equals(result, expected))
        << "result: " << bg::wkt(result) << "\nexpected: " << bg::wkt(expected);
    EXPECT_DOUBLE_EQ(bg::area(result), 84.0);
}
