#include "readData.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

Matrix matrix_fill(int n) {
    Matrix matrix;
    for (int i = 0; i < n; i++) {
        std::vector<double> line;
        for (int j = 0; j < n; j++) {
            line.push_back(std::numeric_limits<double>::infinity());
        }

        matrix.push_back(line);
    }

    return matrix;
}

void get_instance_info(int &dimension, Matrix &cost, std::vector<int> &rnd,
                       const std::string &f_name) {
    std::ifstream f(f_name);
    if (!f.is_open()) {
        std::cerr << "could not open " << f_name << std::endl;
        exit(1);
    }

    std::string line;

    std::getline(f, line);
    dimension = std::stoi(line);
    cost = matrix_fill(dimension);

    for (int i = 0; i < dimension; i++) {
        std::getline(f, line);
        int j = i + 1;
        while (line.find(' ') != std::string::npos) {
            size_t index = line.find(' ');

            cost[i][j] = std::stod(line.substr(0, index));
            cost[j][i] = cost[i][j];

            j += 1;
            line = line.substr(index + 1);
        }
    }

    std::getline(f, line);
    std::getline(f, line);
    std::getline(f, line);
    int rnd_size = std::stoi(line);
    rnd.clear();
    for (int i = 0; i < rnd_size; i++) {
        std::getline(f, line);
        int rnd_value = std::stoi(line);
        rnd.push_back(rnd_value);
    }
}
