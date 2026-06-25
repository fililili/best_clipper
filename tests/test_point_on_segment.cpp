#include "snap_rounding_helper.hpp"
#include <gtest/gtest.h>

using namespace best_clipper;

// Hot pixel Pi(p) = (x-0.5, x+0.5] x (y-0.5, y+0.5]
// Right/top closed, left/bottom open.
// Tests use pixel at (0,0) unless otherwise noted.

// ---------------------------------------------------------------------------
TEST(PointOnSegment, Inside) {
  // through pixel center along axes and diagonals
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
  // through pixel across all 4 corners (from one corner to opposite)
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{2, 2}, point{-2, -2}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-2, 2}, point{2, -2}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-2, -2}, point{2, 2}}));
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{2, -2}, point{-2, 2}}));
}

// ---------------------------------------------------------------------------
TEST(PointOnSegment, Outside) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{2, 2}, point{3, 3}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-3, 0}, point{-2, 0}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{2, -1}, point{2, 1}}));
}

// ---------------------------------------------------------------------------
// Corner grazing: segment touches exactly one pixel corner, then goes away.
// The touch point is in the pixel only when the corner's boundaries are closed
// (right/top: closed; left/bottom: open).
// ---------------------------------------------------------------------------
TEST(PointOnSegment, GrazeRightTopCorner) {
  // (0.5,0.5): right closed, top closed -> included
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{0, 1}, point{1, 0}}));
}

TEST(PointOnSegment, GrazeRightBottomCorner) {
  // (0.5,-0.5): right closed, bottom open -> excluded
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{0, -1}, point{1, 0}}));
}

TEST(PointOnSegment, GrazeLeftTopCorner) {
  // (-0.5,0.5): left open, top closed -> excluded
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 0}, point{0, 1}}));
}

TEST(PointOnSegment, GrazeLeftBottomCorner) {
  // (-0.5,-0.5): left open, bottom open -> excluded
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{0, -1}, point{-1, 0}}));
}

// ---------------------------------------------------------------------------
// Boundary crossing: segment crosses a pixel edge (not at a corner).
// ---------------------------------------------------------------------------
TEST(PointOnSegment, CrossEdges) {
  // bottom edge
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-2, -1}, point{2, 0}}));
  // top edge
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 1}, point{1, 0}}));
  // left edge
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 0}, point{1, 1}}));
  // right edge
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{2, -1}, point{-1, 1}}));
}

// ---------------------------------------------------------------------------
// Edge grazing: segment parallel to a pixel edge, just outside.
// ---------------------------------------------------------------------------
TEST(PointOnSegment, GrazeEdges) {
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, -1}, point{2, -1}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, 1}, point{2, 1}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{-1, -1}, point{-1, 2}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{1, -1}, point{1, 2}}));
}

// ---------------------------------------------------------------------------
TEST(PointOnSegment, ZeroLength) {
  EXPECT_TRUE(
      is_point_on_segment(point{0, 0}, segment{point{0, 0}, point{0, 0}}));
  EXPECT_FALSE(
      is_point_on_segment(point{0, 0}, segment{point{1, 1}, point{1, 1}}));
}

// ---------------------------------------------------------------------------
TEST(PointOnSegment, LargeCoordinates) {
  coordinate_type big = 1000000;
  EXPECT_TRUE(is_point_on_segment(
      point{big, big}, segment{point{big - 1, big}, point{big + 1, big}}));
  EXPECT_FALSE(is_point_on_segment(
      point{big, big}, segment{point{big + 2, big}, point{big + 3, big}}));
}
