#include "include/kdbush.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(-10000, 10000);

    std::vector<int> coords;
    for (std::size_t i = 0; i < 1000000; i++) {
        coords.push_back(dis(gen));
        coords.push_back(dis(gen));
    }

    const auto started = std::chrono::high_resolution_clock::now();
    KDBush<int> index(coords);
    const auto finished = std::chrono::high_resolution_clock::now();

    const auto duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(finished - started).count();

    std::cerr << "indexed 1M points in " << (duration / 1e6) << "ms\n";

    return 0;
}
