#include "Data.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace MLP {

    static std::vector<std::string> ReadAllLines(const std::string & path) {
        std::vector<std::string> lines;
        std::ifstream in(path);
        if (!in) {
            std::cerr << "could not open " << path << std::endl;
            std::exit(1);
        }
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            lines.push_back(line);
        }
        return lines;
    }

    Data::Data() {
        dimension = 0;
        rnd_index = 0;
    }

    void Data::loadData() {
        std::vector<std::string> file = ReadAllLines("../distance_matrix");
        int file_line = 0;
        dimension = std::stoi(file[file_line++]);
        matrix.assign((size_t)dimension * dimension, 0.0);

        for (int i = 1; i < dimension; i++) {
            file_line++;
            int j = i;
            while (file[i].find(" ") != std::string::npos) {
                size_t index = file[i].find(" ");
                matrix[to_1D_matrix(i - 1, j)] = std::stod(file[i].substr(0, index));
                matrix[to_1D_matrix(j, i - 1)] = matrix[to_1D_matrix(i - 1, j)];
                matrix[to_1D_matrix(i - 1, i - 1)] = 0.0;
                file[i] = file[i].substr(index + 1);
                j++;
            }
        }
        matrix[to_1D_matrix(dimension - 1, dimension - 1)] = 0.0;

        file_line++;
        file_line++;
        file_line++;
        int rnd_size = std::stoi(file[file_line]);
        file_line++;
        rnd.assign(rnd_size, 0);
        for (int i = 0; i < rnd_size; i++) {
            rnd[i] = std::stoi(file[file_line++]);
        }

        c.assign((size_t)dimension * dimension, 0.0);
        for (int i = 0; i < dimension; i++) {
            for (int j = i; j < dimension; j++) {
                c[to_1D(i, j)] = getDistance(i, j);
                c[to_1D(j, i)] = getDistance(j, i);
            }
        }
    }

}
