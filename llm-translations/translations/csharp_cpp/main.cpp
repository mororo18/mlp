#include <chrono>
#include <ctime>
#include <cstring>
#include <iostream>

#include "GILS_RVND.hpp"
#include "dbl.hpp"

using namespace MLP;

int main(int argc, char ** argv) {
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
            break;
        }
    }

    if (verbose) {
        std::cout << "Hello World!" << std::endl;
    }
    GILS_RVND tsp;

    std::clock_t cpuStart = std::clock();

    auto start = std::chrono::system_clock::now();
    tsp.solve(verbose);
    auto end = std::chrono::system_clock::now();

    double cpuTime = (double)(std::clock() - cpuStart) / CLOCKS_PER_SEC;

    std::chrono::duration<double> ts = end - start;
    std::cout << "TIME: " << DblToStr(cpuTime) << std::endl;

    std::cout << "wall clock (s): " << DblToStr(ts.count()) << std::endl;

    return 0;
}
