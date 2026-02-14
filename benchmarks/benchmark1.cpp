#include "basic.hpp"

auto benchmark(int size) {
    using namespace std::chrono_literals;
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(-600, 600);

    ring r1, r2;
    point p{ 0, 0 };
    r1.emplace_back(0, 0);
    for (int i = 0; i < size; i++) {
        bg::set<0>(p, bg::get<0>(p) + distrib(gen));
        bg::set<1>(p, bg::get<1>(p) + distrib(gen));
        r1.emplace_back(p);
    }
    r1.emplace_back(0, 0);
    p = point{ 0, 0 };
    r2.emplace_back(0, 0);
    for (int i = 0; i < size; i++) {
        bg::set<0>(p, bg::get<0>(p) + distrib(gen));
        bg::set<1>(p, bg::get<1>(p) + distrib(gen));
        r2.emplace_back(p);
    }
    r2.emplace_back(0, 0);
    std::cout << "\n\n-----------------" << std::endl;
    auto before = std::chrono::system_clock::now();
    std::cout << "run self r1 ----------------------:" << std::endl;
    std::cout << bg::wkt(r1) << std::endl;
    auto sr1 = self_or(r1);
    std::cout << bg::wkt(sr1) << std::endl;
    std::cout << "run self r2 ----------------------:" << std::endl;
    std::cout << bg::wkt(r2) << std::endl;
    auto sr2 = self_or(r2);
    std::cout << bg::wkt(sr2) << std::endl;
    std::cout << "run r1 + r2 ----------------------:" << std::endl;
    auto sum = add(sr1, sr2);
    std::cout << bg::wkt(sum) << std::endl;
    auto after = std::chrono::system_clock::now();
    std::cout << "benchmark size = " << size << ", total runtime: " << (after - before) / 1s << "s" << std::endl;
}

int main()
{
    //while(1)
    benchmark(400);
    /*
    benchmark(100);
    benchmark(200);
    benchmark(400);
    benchmark(1000);
    benchmark(2000);
    benchmark(4000);
    benchmark(10000);
    benchmark(20000);
    benchmark(40000);
    benchmark(80000);
    benchmark(100000);
    benchmark(200000);
    benchmark(400000);
    benchmark(1000000);
    benchmark(2000000);
    benchmark(4000000);
    benchmark(10000000);
    benchmark(20000000);
    benchmark(40000000);
    benchmark(100000000);
    */
}