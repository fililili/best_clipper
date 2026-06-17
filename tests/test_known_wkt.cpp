#include "best_clipper.hpp"
#include <gtest/gtest.h>

using namespace best_clipper;

// ============================================================================
// Intersection tests — verify against known WKT
// ============================================================================

void test_intersection(const std::string &a_wkt, const std::string &b_wkt,
                       const std::string &expected_wkt) {
  multi_polygon a, b, expected;
  bg::read_wkt(a_wkt, a);
  ASSERT_TRUE(bg::is_valid(a));
  bg::read_wkt(b_wkt, b);
  ASSERT_TRUE(bg::is_valid(b));
  bg::read_wkt(expected_wkt, expected);
  ASSERT_TRUE(bg::is_valid(expected));
  auto result = intersection(a, b);
  EXPECT_TRUE(bg::is_valid(result));
  EXPECT_TRUE(bg::equals(result, expected))
      << "Result:   " << bg::wkt(result) << "\n"
      << "Expected: " << bg::wkt(expected);
}

void test_union(const std::string &a_wkt, const std::string &b_wkt,
                const std::string &expected_wkt) {
  multi_polygon a, b, expected;
  bg::read_wkt(a_wkt, a);
  ASSERT_TRUE(bg::is_valid(a));
  bg::read_wkt(b_wkt, b);
  ASSERT_TRUE(bg::is_valid(b));
  bg::read_wkt(expected_wkt, expected);
  ASSERT_TRUE(bg::is_valid(expected));
  auto result = union_(a, b);
  EXPECT_TRUE(bg::is_valid(result));
  EXPECT_TRUE(bg::equals(result, expected))
      << "Result:   " << bg::wkt(result) << "\n"
      << "Expected: " << bg::wkt(expected);
}

void test_robust_self_or(const std::string &a_wkt, const std::string &expected_wkt) {
  multi_polygon a, expected;
  bg::read_wkt(a_wkt, a);
  bg::read_wkt(expected_wkt, expected);
  ASSERT_TRUE(bg::is_valid(expected));
  auto result = robust_self_or(a);
  EXPECT_TRUE(bg::is_valid(result));
  EXPECT_TRUE(bg::equals(result, expected))
      << "Result:   " << bg::wkt(result) << "\n"
      << "Expected: " << bg::wkt(expected);
}

TEST(CoreTest, Intersection) {
  test_intersection("MULTIPOLYGON(((0 0, 0 2, 2 2, 2 0, 0 0)))",
                    "MULTIPOLYGON(((1 1, 1 3, 3 3, 3 1, 1 1)))",
                    "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))");
  test_intersection("MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0)))",
                    "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))",
                    "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))");
  test_intersection("MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))",
                    "MULTIPOLYGON(((2 2, 2 3, 3 3, 3 2, 2 2)))",
                    "MULTIPOLYGON()");
  test_intersection("MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))",
                    "MULTIPOLYGON(((1 0, 1 1, 2 1, 2 0, 1 0)))",
                    "MULTIPOLYGON()");
  test_intersection("MULTIPOLYGON(((0 0, 0 1, 1 1, 1 0, 0 0)))",
                    "MULTIPOLYGON(((1 1, 1 2, 2 2, 2 1, 1 1)))",
                    "MULTIPOLYGON()");
  test_intersection("MULTIPOLYGON(((0 0, 0 10, 10 10, 10 0, 0 0)), ((20 0, 20 "
                    "10, 30 10, 30 0, 20 0)))",
                    "MULTIPOLYGON(((5 5, 5 15, 15 15, 15 5, 5 5)))",
                    "MULTIPOLYGON(((5 5, 5 10, 10 10, 10 5, 5 5)))");
}

TEST(CoreTest, DiagonalSnapRounding) {
  test_union("MULTIPOLYGON(((0 0,3 4,3 0,0 0)))",
             "MULTIPOLYGON(((-1 1,-1 2,0 2,0 1,-1 1)))",
             "MULTIPOLYGON(((0 1,-1 1,-1 2,0 2,0 1)),((0 0,0 1,3 4,3 0,0 0)))");
  test_intersection("MULTIPOLYGON(((0 0,3 4,3 0,0 0)))",
                    "MULTIPOLYGON(((-1 1,-1 2,0 2,0 1,-1 1)))",
                    "MULTIPOLYGON()");

  test_union("MULTIPOLYGON(((0 5,5 10,10 5,5 0,0 5)))",
             "MULTIPOLYGON(((5 5,10 10,15 5,10 0,5 5)))",
             "MULTIPOLYGON(((7 2,5 0,0 5,5 10,7 7,10 10,15 5,10 0,7 2)))");
  test_intersection("MULTIPOLYGON(((0 5,5 10,10 5,5 0,0 5)))",
                    "MULTIPOLYGON(((5 5,10 10,15 5,10 0,5 5)))",
                    "MULTIPOLYGON(((7 2,5 5,7 7,10 5,7 2)))");

  test_union("MULTIPOLYGON(((0 0,0 8,6 0,0 0)))",
             "MULTIPOLYGON(((2 2,2 10,10 10,10 2,2 2)))",
             "MULTIPOLYGON(((6 0,0 0,0 8,2 5,2 10,10 10,10 2,4 2,6 0)))");
  test_intersection("MULTIPOLYGON(((0 0,0 8,6 0,0 0)))",
                    "MULTIPOLYGON(((2 2,2 10,10 10,10 2,2 2)))",
                    "MULTIPOLYGON(((4 2,2 2,2 5,4 2)))");

  test_union("MULTIPOLYGON(((0 0,0 10,10 10,10 0,0 0)))",
             "MULTIPOLYGON(((3 9,6 11,13 7,10 5,3 9)))",
             "MULTIPOLYGON(((10 5,10 0,0 0,0 10,4 10,6 11,8 10,10 10,10 9,13 "
             "7,10 5)))");
  test_intersection("MULTIPOLYGON(((0 0,0 10,10 10,10 0,0 0)))",
                    "MULTIPOLYGON(((3 9,6 11,13 7,10 5,3 9)))",
                    "MULTIPOLYGON(((10 5,3 9,4 10,8 10,10 9,10 5)))");

  test_union("MULTIPOLYGON(((0 0,0 4,8 4,8 0,0 0)))",
             "MULTIPOLYGON(((4 2,4 6,12 6,12 2,4 2)))",
             "MULTIPOLYGON(((8 2,8 0,0 0,0 4,4 4,4 6,12 6,12 2,8 2)))");
  test_intersection("MULTIPOLYGON(((0 0,0 4,8 4,8 0,0 0)))",
                    "MULTIPOLYGON(((4 2,4 6,12 6,12 2,4 2)))",
                    "MULTIPOLYGON(((8 2,4 2,4 4,8 4,8 2)))");

  test_robust_self_or("MULTIPOLYGON(((0 3,3 6,6 3,3 0,0 3)),"
               "((3 3,6 6,9 3,6 0,3 3)))",
               "MULTIPOLYGON(((4 1,3 0,0 3,3 6,4 4,6 6,9 3,6 0,4 1)))");
}

TEST(CoreTest, Union) {
  test_union("MULTIPOLYGON(((-59 867,-36 492,-182 486,-59 867)))",
             "MULTIPOLYGON(((-220 877,-54 821,-402 541,-808 638,-220 877)))",
             "MULTIPOLYGON(((-220 877,-72 827,-59 867,-56 822,-54 821,-56 "
             "819,-36 492,-182 486,-81 799,-402 541,-808 638,-220 877)))");
  test_union(
      "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))",
      "MULTIPOLYGON(((2 2, 2 4, 4 4, 4 2, 2 2)))",
      "MULTIPOLYGON(((0 0, 0 3, 2 3, 2 4, 4 4, 4 2, 3 2, 3 0, 0 0), (1 1, 2 1, "
      "2 2, 1 2, 1 1)))");
  test_union("MULTIPOLYGON(((0 0, 0 2, 5 1, 5 0, 0 0)))",
             "MULTIPOLYGON(((3 1, 0 2, 5 1, 3 1)))",
             "MULTIPOLYGON(((0 2,3 1,5 1,5 0,0 0,0 2)))");
  test_union("MULTIPOLYGON(((-1 -1, -1 3, 3 3, 3 -1, -1 -1), (0 0, 2 0, 2 2, 0 "
             "2, 0 0)))",
             "MULTIPOLYGON(((1 1, 1 4, 4 4, 4 1, 1 1)))",
             "MULTIPOLYGON(((-1 3,1 3,1 4,4 4,4 1,3 1,3 -1,-1 -1,-1 3),(2 0,2 "
             "1,1 1,1 2,0 2,0 0,2 0)))");
  test_union("MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 "
             "1), (6 6, 8 6, 8 8, 6 8, 6 6)))",
             "MULTIPOLYGON(((2 2, 2 7, 7 7, 7 2, 2 2)))",
             "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 2, 2 2, 2 "
             "3, 1 3, 1 1), (8 8, 6 8, 6 7, 7 7, 7 6, 8 6, 8 8)))");
  test_union(
      "MULTIPOLYGON(((0 0, 1 1, 2 1, 2 2, 3 3, 3 0, 0 0)))",
      "MULTIPOLYGON(((0 0, 0 3, 3 3, 2 2, 1 2, 1 1, 0 0)))",
      "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))");
  test_union(
      "MULTIPOLYGON(((-1461 -786,-1417 -833,-1389 -830,-1450 -775,-1061 "
      "-372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 "
      "-741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 "
      "-773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 "
      "-1136,-956 -1205,-1030 -1246,-1608 -939,-1461 -786),(-1058 -1133,-1244 "
      "-963,-1301 -1039,-1058 -1133)),((-578 -670,-688 -679,-802 -443,-826 "
      "-400,-578 -670)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 "
      "-666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 "
      "-1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 "
      "-838,-378 -839)))",
      "MULTIPOLYGON(((-1450 -1280, -1450 -800, -1200 -1000, -1000 -1280, -1450 "
      "-1280)))",
      "MULTIPOLYGON(((-1461 -786,-1442 -807,-1410 -832,-1389 -830,-1450 "
      "-775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 "
      "-855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 "
      "-790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 "
      "-1136,-956 -1205,-1026 -1244,-1000 -1280,-1450 -1280,-1450 -1023,-1608 "
      "-939,-1461 -786),(-1245 -964,-1228 -977,-1244 -963,-1245 -964),(-1123 "
      "-1108,-1058 -1133,-1193 -1009,-1123 -1108)),((-826 -400,-578 -670,-688 "
      "-679,-802 -443,-826 -400)),((-433 -621,-300 -283,-289 -294,-432 "
      "-660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 "
      "-1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 "
      "-844,-380 -838,-378 -839)))");
}

TEST(CoreTest, SnapRoundingRegression) {
  // These two overlapping 11-gons trigger a snap-rounding edge case that
  // caused power-balance assertion failures before commit 4a3c49d.
  // The expected union result is a multi-polygon whose components touch;
  // we compare directly without bg::is_valid on the expected value.
  const char *a_wkt =
      "MULTIPOLYGON(((633 255,626 248,619 245,609 251,603 257,601 268,611 "
      "276,616 283,626 281,633 274,636 265,633 255)))";
  const char *b_wkt =
      "MULTIPOLYGON(((625 263,622 259,617 261,612 262,611 267,611 269,616 "
      "275,616 275,624 275,627 271,626 267,625 263)))";
  const char *u_wkt =
      "MULTIPOLYGON(((636 265,633 255,626 248,619 245,609 251,603 257,601 "
      "268,611 276,616 283,626 281,633 274,636 265)),((611 267,611 269,616 "
      "275,624 275,627 271,626 267,625 263,622 259,617 261,612 262,611 267)))";
  const char *i_wkt = "MULTIPOLYGON()";

  multi_polygon a, b, u_exp, i_exp, u_res, i_res;
  bg::read_wkt(a_wkt, a);
  ASSERT_TRUE(bg::is_valid(a));
  bg::read_wkt(b_wkt, b);
  ASSERT_TRUE(bg::is_valid(b));
  bg::read_wkt(u_wkt, u_exp);
  bg::read_wkt(i_wkt, i_exp);

  u_res = union_(a, b);
  EXPECT_TRUE(bg::equals(u_res, u_exp))
      << "Union result:   " << bg::wkt(u_res) << "\n"
      << "Expected: " << bg::wkt(u_exp);

  i_res = intersection(a, b);
  EXPECT_TRUE(bg::equals(i_res, i_exp))
      << "Intersection result:   " << bg::wkt(i_res) << "\n"
      << "Expected: " << bg::wkt(i_exp);
}
