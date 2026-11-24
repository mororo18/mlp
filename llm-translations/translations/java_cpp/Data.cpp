#include "Data.h"

#include <fstream>
#include <iostream>
#include <string>

Data::Data() {
    dimension = 0;
    rnd_size = 0;
}

void Data::loadData() {
    std::string file_name = "../distance_matrix";

    std::ifstream f_reader(file_name);
    if (!f_reader.is_open()) {
        std::cout << "Error while reading " << file_name << " file" << std::endl;
        return;
    }

    std::string line;
    std::getline(f_reader, line);

    dimension = std::stoi(line);
    matrix.assign(dimension, std::vector<double>(dimension, 0.0));
    for (int k = 0; k < dimension; k++) matrix[k][k] = 0.0;

    for (int i = 0; i < dimension; i++) {
        std::getline(f_reader, line);

        int j = i + 1;
        while (line.find(' ') != std::string::npos) {
            std::string cost = line.substr(0, line.find(' '));
            line = line.substr(line.find(' ') + 1, line.length());

            double cost_value = std::stod(cost);
            matrix[i][j] = cost_value;
            matrix[j][i] = matrix[i][j];
            j++;
        }
    }

    std::getline(f_reader, line);
    std::getline(f_reader, line);
    std::getline(f_reader, line);

    rnd_size = std::stoi(line);
    rnd.assign(rnd_size, 0);

    for (int i = 0; i < rnd_size; i++) {
        std::getline(f_reader, line);
        rnd[i] = std::stoi(line);
    }

    f_reader.close();
}
