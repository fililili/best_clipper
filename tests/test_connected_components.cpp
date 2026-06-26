#include "best_clipper.hpp"
#include "graph_helper.hpp"
#include <gtest/gtest.h>

using namespace best_clipper;

// ============================================================================
// connected_components unit tests
// ============================================================================

TEST(ConnectedComponents, EmptyGraph) {
  auto [comp, sizes] = connected_components(5, {});
  EXPECT_EQ(comp.size(), 5u);
  EXPECT_EQ(sizes.size(), 5u);
  for (std::size_t i = 0; i < 5; i++) {
    EXPECT_EQ(comp[i], i) << "vertex " << i;
    EXPECT_EQ(sizes[comp[i]], 1u);
  }
}

TEST(ConnectedComponents, SingleEdge) {
  auto [comp, sizes] = connected_components(4, {{0, 1}});
  EXPECT_EQ(comp[0], comp[1]);
  EXPECT_NE(comp[0], comp[2]);
  EXPECT_NE(comp[0], comp[3]);
  EXPECT_NE(comp[2], comp[3]);
  EXPECT_EQ(sizes.size(), 3u);
  EXPECT_EQ(sizes[comp[0]], 2u);
  EXPECT_EQ(sizes[comp[2]], 1u);
  EXPECT_EQ(sizes[comp[3]], 1u);
}

TEST(ConnectedComponents, Chain) {
  auto [comp, sizes] = connected_components(5, {{0, 1}, {1, 2}, {2, 3}});
  EXPECT_EQ(comp[0], comp[1]);
  EXPECT_EQ(comp[0], comp[2]);
  EXPECT_EQ(comp[0], comp[3]);
  EXPECT_NE(comp[0], comp[4]);
  EXPECT_EQ(sizes.size(), 2u);
  EXPECT_EQ(sizes[comp[0]], 4u);
  EXPECT_EQ(sizes[comp[4]], 1u);
}

TEST(ConnectedComponents, Disconnected) {
  auto [comp, sizes] = connected_components(6, {{0, 1}, {2, 3}, {4, 5}});
  EXPECT_EQ(comp[0], comp[1]);
  EXPECT_EQ(comp[2], comp[3]);
  EXPECT_EQ(comp[4], comp[5]);
  EXPECT_NE(comp[0], comp[2]);
  EXPECT_NE(comp[0], comp[4]);
  EXPECT_NE(comp[2], comp[4]);
  EXPECT_EQ(sizes.size(), 3u);
  for (auto s : sizes)
    EXPECT_EQ(s, 2u);
}

TEST(ConnectedComponents, DuplicateEdges) {
  auto [comp, sizes] = connected_components(3, {{0, 1}, {0, 1}, {1, 0}});
  EXPECT_EQ(comp[0], comp[1]);
  EXPECT_NE(comp[0], comp[2]);
  EXPECT_EQ(sizes.size(), 2u);
  EXPECT_EQ(sizes[comp[0]], 2u);
  EXPECT_EQ(sizes[comp[2]], 1u);
}

TEST(ConnectedComponents, SelfLoop) {
  auto [comp, sizes] = connected_components(3, {{0, 0}, {0, 1}});
  EXPECT_EQ(comp[0], comp[1]);
  EXPECT_NE(comp[0], comp[2]);
  EXPECT_EQ(sizes.size(), 2u);
  EXPECT_EQ(sizes[comp[0]], 2u);
  EXPECT_EQ(sizes[comp[2]], 1u);
}

TEST(ConnectedComponents, StarGraph) {
  auto [comp, sizes] =
      connected_components(5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}});
  for (std::size_t i = 1; i < 5; i++)
    EXPECT_EQ(comp[0], comp[i]);
  EXPECT_EQ(sizes.size(), 1u);
  EXPECT_EQ(sizes[0], 5u);
}

TEST(ConnectedComponents, IsolatedVertex) {
  auto [comp, sizes] = connected_components(4, {{0, 1}});
  EXPECT_EQ(comp[0], comp[1]);
  EXPECT_NE(comp[2], comp[0]);
  EXPECT_NE(comp[3], comp[0]);
  EXPECT_NE(comp[2], comp[3]);
  EXPECT_EQ(sizes.size(), 3u);
  EXPECT_EQ(sizes[comp[0]], 2u);
}

TEST(ConnectedComponents, LargeChain) {
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  for (std::size_t i = 0; i + 1 < 1000; i++)
    edges.emplace_back(i, i + 1);
  auto [comp, sizes] = connected_components(1000, edges);
  for (std::size_t i = 1; i < 1000; i++)
    EXPECT_EQ(comp[0], comp[i]);
  EXPECT_EQ(comp[0], 0u);
  EXPECT_EQ(sizes.size(), 1u);
  EXPECT_EQ(sizes[0], 1000u);
}
