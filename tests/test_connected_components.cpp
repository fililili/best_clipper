#include <gtest/gtest.h>
#include "uint32_adaptor.hpp"

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
