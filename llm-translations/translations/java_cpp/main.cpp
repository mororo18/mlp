#include "GILS_RVND.h"

#include <chrono>
#include <ctime>
#include <iostream>
#include <string>

static long processCpuTime() {
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return (long) ts.tv_sec * 1000000000L + (long) ts.tv_nsec;
}

static long nanoTime() {
    return (long) std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void read(bool verbose) {
    GILS_RVND tsp = GILS_RVND(verbose);
    long cpuTime = - processCpuTime();
    long time = - nanoTime();
    tsp.solve();
    time += nanoTime();
    cpuTime += processCpuTime();

    std::cout << "TIME: ";
    std::cout << cpuTime/10e8 << std::endl;
    std::cout << "wall clock (s): ";
    std::cout << time/10e8 << std::endl;
}

int main(int argc, char *argv[]) {
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
    }
    read(verbose);
    return 0;
}
