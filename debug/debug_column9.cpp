#include <cstdio>
#include "core.hpp"
using namespace best_clipper;

ring make_rect(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2) {
    ring r;
    r.push_back(point{x1, y1});
    r.push_back(point{x1, y2});
    r.push_back(point{x2, y2});
    r.push_back(point{x2, y1});
    r.push_back(point{x1, y1});
    return r;
}

ring make_hole(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2) {
    ring r;
    r.push_back(point{x1, y1});
    r.push_back(point{x2, y1});
    r.push_back(point{x2, y2});
    r.push_back(point{x1, y2});
    r.push_back(point{x1, y1});
    return r;
}

polygon make_rect_poly(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2) {
    polygon p;
    p.outer() = make_rect(x1, y1, x2, y2);
    return p;
}

int main() {
    multi_polygon a;
    {
        polygon p; p.outer() = make_rect(0, 0, 100, 100);
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 10; j++)
                p.inners().push_back(make_hole(2 + i * 10, 2 + j * 10, 8 + i * 10, 8 + j * 10));
        for (int j = 0; j < 10; j++)
            p.inners().push_back(make_hole(92, 2 + j * 10, 98, 8 + j * 10));
        a.push_back(p);
    }
    multi_polygon b;
    b.push_back(make_rect_poly(1, 1, 52, 99));

    std::printf("Input: %zu holes total\n",
                a[0].inners().size());

    auto result = add(a, b);
    std::printf("Result area: %.0f\n", bg::area(result));
    std::printf("Result valid: %s\n", bg::is_valid(result) ? "yes" : "no");
    std::printf("Result polys: %zu\n", result.size());
    return 0;
}
