#include <cstdio>
#include <sstream>
#include "core.hpp"
using namespace best_clipper;

int main() {
    // Check area conventions
    // Ring going (0,0)->(0,2)->(2,2)->(2,0)->(0,0): traditional CCW
    ring r1;
    bg::read_wkt("((0 0,0 2,2 2,2 0,0 0))", r1);
    printf("Ring1 area: %.0f (should be positive for CCW)\n", bg::area(r1));

    // Ring going (0,0)->(2,0)->(2,2)->(0,2)->(0,0): traditional CW
    ring r2;
    bg::read_wkt("((0 0,2 0,2 2,0 2,0 0))", r2);
    printf("Ring2 area: %.0f (should be negative for CW)\n", bg::area(r2));

    // What our algorithm produces for single square
    ring r3;
    bg::read_wkt("((2 2,0 2,0 0,2 0,2 2))", r3);
    printf("Ring3 area: %.0f\n", bg::area(r3));

    // Check polygon convention: CCW outer gives positive area
    polygon p_ccw;
    bg::read_wkt("((0 0,0 2,2 2,2 0,0 0))", p_ccw);
    printf("Polygon (0,0)-(0,2)-(2,2)-(2,0) area: %.0f, valid: %d\n",
           bg::area(p_ccw), bg::is_valid(p_ccw));

    // Now test our algorithm with poly explicitly CW
    multi_polygon a;
    bg::read_wkt("MULTIPOLYGON(((0 0,2 0,2 2,0 2,0 0)))", a);
    printf("CW poly area: %.0f, valid: %d\n", bg::area(a), bg::is_valid(a));
    auto result = self_or(a);
    std::stringstream ss;
    ss << bg::wkt(result);
    printf("Result: %s\n", ss.str().c_str());
    printf("Result area: %.0f, valid: %d\n", bg::area(result), bg::is_valid(result));

    return 0;
}
