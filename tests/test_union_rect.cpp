#include "best_clipper.hpp"
#include <gtest/gtest.h>

using namespace best_clipper;

// ============================================================================
// Rectangle union / self_or — many overlapping squares
// ============================================================================

void test_union_rectangle(int size) {
  multi_polygon a, b;
  for (int i = 0; i < size; i++) {
    a.emplace_back(polygon{{{0 + 2 * i, 0 + 2 * i},
                            {0 + 2 * i, 2 + 2 * i},
                            {2 + 2 * i, 2 + 2 * i},
                            {2 + 2 * i, 0 + 2 * i},
                            {0 + 2 * i, 0 + 2 * i}}});
    b.emplace_back(polygon{{{1 + 2 * i, 1 + 2 * i},
                            {1 + 2 * i, 3 + 2 * i},
                            {3 + 2 * i, 3 + 2 * i},
                            {3 + 2 * i, 1 + 2 * i},
                            {1 + 2 * i, 1 + 2 * i}}});
  }
  EXPECT_TRUE(bg::is_valid(a));
  EXPECT_TRUE(bg::is_valid(b));
  EXPECT_TRUE(bg::equals(self_or(a), a));
  EXPECT_TRUE(bg::equals(self_or(b), b));
  auto result = union_(a, b);
  EXPECT_TRUE(bg::is_valid(result));
  EXPECT_DOUBLE_EQ(bg::area(result), 1 + 6 * size);
}

void test_self_or_rectangle(int size) {
  multi_polygon poly;
  for (int i = 0; i < size; i++) {
    poly.emplace_back(polygon{{{0 + i, 0 + i},
                               {0 + i, 2 + i},
                               {2 + i, 2 + i},
                               {2 + i, 0 + i},
                               {0 + i, 0 + i}}});
  }
  auto result = self_or(poly);
  EXPECT_TRUE(bg::is_valid(result));
  EXPECT_DOUBLE_EQ(bg::area(result), 1 + 3 * size);
}

TEST(CoreUnion, RectangleUnion) {
  test_union_rectangle(100);
  test_union_rectangle(3);
  test_union_rectangle(1010);
  test_union_rectangle(300);
  test_union_rectangle(2521);
}

TEST(CoreUnion, RectangleSelfOr) {
  test_self_or_rectangle(100);
  test_self_or_rectangle(3);
  test_self_or_rectangle(1010);
  test_self_or_rectangle(300);
  test_self_or_rectangle(2521);
}
