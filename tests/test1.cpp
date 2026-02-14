#include "basic.hpp"

void test_union(std::string first_s, std::string second_s, std::string ret_s) {
    multi_polygon first, second, ret;
    bg::read_wkt(first_s, first);
    assert(bg::is_valid(first));
    bg::read_wkt(second_s, second);
    assert(bg::is_valid(second));
    bg::read_wkt(ret_s, ret);
    assert(bg::is_valid(ret));
    std::cout << bg::wkt(add(first, second)) << std::endl;
    assert(bg::equals(add(first, second), ret));
}

void test_union_rectangle(int size) {
    using namespace std::chrono_literals;
    multi_polygon first, second;
    for (int i = 0; i < size; i++) {
        first.emplace_back(polygon{ {{0 + 2 * i, 0 + 2 * i}, {0 + 2 * i, 2 + 2 * i}, {2 + 2 * i, 2 + 2 * i}, {2 + 2 * i, 0 + 2 * i}, {0 + 2 * i, 0 + 2 * i}} });
        second.emplace_back(polygon{ {{1 + 2 * i, 1 + 2 * i}, {1 + 2 * i, 3 + 2 * i}, {3 + 2 * i, 3 + 2 * i}, {3 + 2 * i, 1 + 2 * i}, {1 + 2 * i, 1 + 2 * i}} });
    }
    auto before = std::chrono::system_clock::now();
    assert(bg::is_valid(first));
    assert(bg::is_valid(second));
    assert(bg::equals(self_or(first), first));
    assert(bg::equals(self_or(second), second));
    auto ret = add(first, second);
    assert(bg::is_valid(ret));
    assert(bg::area(ret) == (1 + 6 * size));
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

void test_self_or_rectangle(int size) {
    using namespace std::chrono_literals;
    multi_polygon poly;
    for (int i = 0; i < size; i++) {
        poly.emplace_back(polygon{ {{0 + i, 0 + i}, {0 + i, 2 + i}, {2 + i, 2 + i}, {2 + i, 0 + i}, {0 + i, 0 + i}} });
    }
    auto before = std::chrono::system_clock::now();
    std::cout << bg::wkt(poly) << std::endl;
    auto ret = self_or(poly);
    assert(bg::is_valid(ret));
    assert(bg::area(ret) == (1 + 3 * size));
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

int main()
{
    multi_polygon one, two, ret;

    test_union(
        "MULTIPOLYGON(((-59 867,-36 492,-182 486,-59 867)))",
        "MULTIPOLYGON(((-220 877,-54 821,-402 541,-808 638,-220 877)))",
        "MULTIPOLYGON(((-220 877,-72 827,-59 867,-56 822,-54 821,-56 819,-36 492,-182 486,-81 799,-402 541,-808 638,-220 877)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))",
        "MULTIPOLYGON(((2 2, 2 4, 4 4, 4 2, 2 2)))",
        "MULTIPOLYGON(((0 0, 0 3, 2 3, 2 4, 4 4, 4 2, 3 2, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 2, 5 1, 5 0, 0 0)))",
        "MULTIPOLYGON(((3 1, 0 2, 5 1, 3 1)))",
        "MULTIPOLYGON(((0 2,3 1,5 1,5 0,0 0,0 2)))"
    );
    test_union(
        "MULTIPOLYGON(((-1 -1, -1 3, 3 3, 3 -1, -1 -1), (0 0, 2 0, 2 2, 0 2, 0 0)))",
        "MULTIPOLYGON(((1 1, 1 4, 4 4, 4 1, 1 1)))",
        "MULTIPOLYGON(((-1 3,1 3,1 4,4 4,4 1,3 1,3 -1,-1 -1,-1 3),(2 0,2 1,1 1,1 2,0 2,0 0,2 0)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1), (6 6, 8 6, 8 8, 6 8, 6 6)))",
        "MULTIPOLYGON(((2 2, 2 7, 7 7, 7 2, 2 2)))",
        "MULTIPOLYGON(((0 0, 0 9, 9 9, 9 0, 0 0), (1 1, 3 1, 3 2, 2 2, 2 3, 1 3, 1 1), (8 8, 6 8, 6 7, 7 7, 7 6, 8 6, 8 8)))"
    );
    test_union(
        "MULTIPOLYGON(((0 0, 1 1, 2 1, 2 2, 3 3, 3 0, 0 0)))",
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 2 2, 1 2, 1 1, 0 0)))",
        "MULTIPOLYGON(((0 0, 0 3, 3 3, 3 0, 0 0), (1 1, 2 1, 2 2, 1 2, 1 1)))"
    );
    test_union(
        "MULTIPOLYGON(((-1461 -786,-1417 -833,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1030 -1246,-1608 -939,-1461 -786),(-1058 -1133,-1244 -963,-1301 -1039,-1058 -1133)),((-578 -670,-688 -679,-802 -443,-826 -400,-578 -670)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))",
        "MULTIPOLYGON(((-1450 -1280, -1450 -800, -1200 -1000, -1000 -1280, -1450 -1280)))",
        "MULTIPOLYGON(((-1461 -786,-1442 -807,-1410 -832,-1389 -830,-1450 -775,-1061 -372,-720 -681,-1007 -702,-1005 -642,-1145 -830,-873 -855,-658 -741,-660 -736,-656 -740,-561 -689,-535 -717,-497 -747,-634 -790,-642 -773,-666 -800,-849 -858,-748 -867,-807 -964,-1012 -1200,-913 -1136,-956 -1205,-1026 -1244,-1000 -1280,-1450 -1280,-1450 -1023,-1608 -939,-1461 -786),(-1245 -964,-1228 -977,-1244 -963,-1245 -964),(-1123 -1108,-1058 -1133,-1193 -1009,-1123 -1108)),((-826 -400,-578 -670,-688 -679,-802 -443,-826 -400)),((-433 -621,-300 -283,-289 -294,-432 -660,-517 -666,-433 -621)),((-178 -1224,-344 -907,-361 -853,-84 -1068,273 -1394,61 -1669,-438 -1425,-178 -1224)),((-378 -839,-376 -844,-380 -838,-378 -839)))"
    );
    test_union_rectangle(100);

    return 0;
}
