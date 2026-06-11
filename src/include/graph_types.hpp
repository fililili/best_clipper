#pragma once

#include <cstdint>

namespace best_clipper {
struct edge_t { std::size_t start; std::size_t end; };
struct edge_with_power_t { std::size_t start, end; int power; };
}