#include "snap_rounding_helper.hpp"
#include <gtest/gtest.h>

using namespace best_clipper;

// Hot pixel Pi(p) = (x-0.5, x+0.5] x (y-0.5, y+0.5]
// Right/top closed, left/bottom open.
// Snap rounding in get_intersection() rounds half-integers toward -inf.

// ============================================================================
// is_point_on_segment
// ============================================================================

TEST(SnapRounding, PointOnSegmentInside) {
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 0}, point{1, 0}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{0, -1}, point{0, 1}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, -1}, point{1, 1}}));
  // pixel center as endpoint
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{0, 0}, point{1, 1}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, -1}, point{0, 0}}));
  // through pixel across all 4 corners
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{2, 2}, point{-2, -2}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-2, 2}, point{2, -2}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-2, -2}, point{2, 2}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{2, -2}, point{-2, 2}}));
}

TEST(SnapRounding, PointOnSegmentOutside) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{2, 2}, point{3, 3}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-3, 0}, point{-2, 0}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{2, -1}, point{2, 1}}));
}

// Corner grazing: segment touches exactly one pixel corner, then goes away.
// The touch point is in the pixel only when the corner's boundaries are closed
// (right/top: closed; left/bottom: open).
TEST(SnapRounding, GrazeRightTopCorner) {
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{0, 1}, point{1, 0}}));
}

TEST(SnapRounding, GrazeRightBottomCorner) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{0, -1}, point{1, 0}}));
}

TEST(SnapRounding, GrazeLeftTopCorner) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 0}, point{0, 1}}));
}

TEST(SnapRounding, GrazeLeftBottomCorner) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{0, -1}, point{-1, 0}}));
}

TEST(SnapRounding, CrossEdges) {
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-2, -1}, point{2, 0}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 1}, point{1, 0}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 0}, point{1, 1}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{2, -1}, point{-1, 1}}));
}

TEST(SnapRounding, GrazeEdges) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, -1}, point{2, -1}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 1}, point{2, 1}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, -1}, point{-1, 2}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{1, -1}, point{1, 2}}));
}

TEST(SnapRounding, PointOnSegmentZeroLength) {
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{0, 0}, point{0, 0}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{1, 1}, point{1, 1}}));
}

TEST(SnapRounding, PointOnSegmentLargeCoordinates) {
  coordinate_type big = 1000000;
  EXPECT_TRUE(is_point_on_segment(
      point{big, big}, segment{point{big - 1, big}, point{big + 1, big}}));
  EXPECT_FALSE(is_point_on_segment(
      point{big, big}, segment{point{big + 2, big}, point{big + 3, big}}));
}

// ============================================================================
// get_intersection
// ============================================================================

TEST(SnapRounding, IntersectionBasic) {
  auto r = get_intersection(segment{point{0, 0}, point{2, 2}},
                            segment{point{2, 0}, point{0, 2}});
  ASSERT_TRUE(r.has_value());
  EXPECT_EQ(bg::get<0>(*r), 1);
  EXPECT_EQ(bg::get<1>(*r), 1);
}

TEST(SnapRounding, IntersectionSnapHalfInteger) {
  // true intersection at (1.5, 1.5). Snap half-down gives (1, 1).
  auto r = get_intersection(segment{point{0, 0}, point{3, 3}},
                            segment{point{3, 0}, point{0, 3}});
  ASSERT_TRUE(r.has_value());
  EXPECT_EQ(bg::get<0>(*r), 1);
  EXPECT_EQ(bg::get<1>(*r), 1);
}

TEST(SnapRounding, IntersectionAxisAligned) {
  // vertical crossing horizontal
  auto r = get_intersection(segment{point{2, -1}, point{2, 3}},
                            segment{point{0, 1}, point{4, 1}});
  ASSERT_TRUE(r.has_value());
  EXPECT_EQ(bg::get<0>(*r), 2);
  EXPECT_EQ(bg::get<1>(*r), 1);
}

TEST(SnapRounding, IntersectionNone) {
  // parallel
  EXPECT_FALSE(get_intersection(segment{point{0, 0}, point{2, 2}},
                                segment{point{1, 0}, point{3, 2}})
                   .has_value());
  // don't cross
  EXPECT_FALSE(get_intersection(segment{point{0, 0}, point{1, 1}},
                                segment{point{2, 2}, point{3, 3}})
                   .has_value());
  // endpoint touch
  EXPECT_FALSE(get_intersection(segment{point{0, 0}, point{1, 1}},
                                segment{point{1, 1}, point{2, 0}})
                   .has_value());
}
