#include "Data.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

void get_instance_info(int &dimension, std::vector<double> &matrix,
                       std::vector<int> &rand_values) {

    std::ifstream f("../distance_matrix");
    if (!f.is_open()) {
        std::fprintf(stderr, "could not open ../distance_matrix\n");
        std::exit(1);
    }

    std::string line;
    std::getline(f, line);
    dimension = std::stoi(line);

    // 1-based indexing: index (i, j) lives at i*(dimension+1) + j
    matrix.assign((dimension + 1) * (dimension + 1), 0.0);

    rand_values.clear();
    bool rnd = false;
    int sz = -1;
    bool name = false;
    std::string inst = "";

    int i = 1;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        if (sz > 0) {
            rand_values.push_back(std::stoi(line));
            continue;
        }

        if (rnd == true) {
            sz = std::stoi(line);
            continue;
        }

        if (line == "RND") {
            rnd = true;
            continue;
        }

        if (name == true) {
            inst = line;
            name = false;
            continue;
        }

        if (line == "EOF") {
            name = true;
            continue;
        }

        int j = i + 1;
        std::istringstream ss(line);
        double d;
        while (ss >> d) {
            matrix[i * (dimension + 1) + j] = d;
            matrix[j * (dimension + 1) + i] = d;

            j += 1;
        }

        i += 1;
    }

    (void)inst;

    f.close();
}
