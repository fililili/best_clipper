// Debug: why does union of square-with-hole + square lose the hole?
#include <iostream>
#include "uint32_adaptor.hpp"

int main() {
    multi_polygon_s32 a, b;
    bg::read_wkt("MULTIPOLYGON(((-1 -1, -1 3, 3 3, 3 -1, -1 -1), (0 0, 2 0, 2 2, 0 2, 0 0)))", a);
    bg::read_wkt("MULTIPOLYGON(((1 1, 1 4, 4 4, 4 1, 1 1)))", b);

    std::cout << "a is_valid: " << bg::is_valid(a) << ", area: " << bg::area(a) << std::endl;
    std::cout << "b is_valid: " << bg::is_valid(b) << ", area: " << bg::area(b) << std::endl;

    // Check Boost's ground truth
    multi_polygon_s32 bg_result;
    bg::union_(a, b, bg_result);
    std::cout << "Boost union: " << bg::wkt(bg_result) << std::endl;
    std::cout << "Boost is_valid: " << bg::is_valid(bg_result) << std::endl;

    auto result = add(a, b);
    std::cout << "Our union:   " << bg::wkt(result) << std::endl;
    std::cout << "Our is_valid: " << bg::is_valid(result) << std::endl;
    std::cout << "Our area: " << bg::area(result) << std::endl;
    std::cout << "Boost area: " << bg::area(bg_result) << std::endl;

    // Also test with shifted coords to rule out negative values
    multi_polygon_s32 a2, b2;
    bg::read_wkt("MULTIPOLYGON(((1999 1999, 1999 2003, 2003 2003, 2003 1999, 1999 1999), (2000 2000, 2002 2000, 2002 2002, 2000 2002, 2000 2000)))", a2);
    bg::read_wkt("MULTIPOLYGON(((2001 2001, 2001 2004, 2004 2004, 2004 2001, 2001 2001)))", b2);
    auto result2 = add(a2, b2);
    std::cout << "Shifted union: " << bg::wkt(result2) << std::endl;
    std::cout << "Shifted is_valid: " << bg::is_valid(result2) << std::endl;

    return 0;
}
