#include "core.hpp"
#include <iostream>
using namespace best_clipper;

int main() {
    multi_polygon first, second, expected;
    bg::read_wkt("MULTIPOLYGON(((-59 867,-36 492,-182 486,-59 867)))", first);
    bg::read_wkt("MULTIPOLYGON(((-220 877,-54 821,-402 541,-808 638,-220 877)))", second);
    bg::read_wkt("MULTIPOLYGON(((-220 877,-72 827,-59 867,-56 822,-54 821,-56 819,-36 492,-182 486,-81 799,-402 541,-808 638,-220 877)))", expected);
    
    std::cerr << "First valid: " << bg::is_valid(first) << " area=" << bg::area(first) << std::endl;
    std::cerr << "Second valid: " << bg::is_valid(second) << " area=" << bg::area(second) << std::endl;
    
    auto result = add(first, second);
    std::cerr << "Result valid: " << bg::is_valid(result) << " area=" << bg::area(result) << std::endl;
    std::cerr << "Result WKT: " << bg::wkt(result) << std::endl;
    std::cerr << "Expected WKT: " << bg::wkt(expected) << std::endl;
    std::cerr << "Expected valid: " << bg::is_valid(expected) << " area=" << bg::area(expected) << std::endl;
    std::cerr << "Equals: " << bg::equals(result, expected) << std::endl;
    
    // Check area difference
    std::cerr << "Area diff: " << (bg::area(result) - bg::area(expected)) << std::endl;
    
    return 0;
}
