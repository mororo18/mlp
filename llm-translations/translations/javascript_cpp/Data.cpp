#include "Data.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

static std::vector<std::string> read_file_lines(const char * path) {
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string text = buffer.str();

    std::vector<std::string> file;
    std::string line;
    std::istringstream stream(text);
    while (std::getline(stream, line, '\n')) {
        file.push_back(line);
    }
    return file;
}

static std::vector<std::string> split_spaces(const std::string & line) {
    std::vector<std::string> parts;
    std::string token;
    std::istringstream stream(line);
    while (std::getline(stream, token, ' ')) {
        parts.push_back(token);
    }
    return parts;
}

static void matrix_fill(std::vector<std::vector<double> > & matrix, int n) {
    matrix.assign(n, std::vector<double>(n));
}

InfoLoadResult info_load(std::vector<std::vector<double> > & c) {
    std::vector<std::string> file = read_file_lines("../distance_matrix");

    int index = 0;

    int dimension = std::atoi(file[index++].c_str());
    matrix_fill(c, dimension);

    for (int i = 0; i < dimension; i++) {
        c[i][i] = 0.0;
        std::vector<std::string> cost = split_spaces(file[index]);
        for (int j = i + 1; j < dimension; j++) {
            c[i][j] = std::atof(cost[j - (i + 1)].c_str());
            c[j][i] = std::atof(cost[j - (i + 1)].c_str());
        }

        index++;
    }

    index++;
    index++;
    int rnd_size = std::atoi(file[index++].c_str());
    std::vector<int> rnd(rnd_size);

    for (int i = 0; i < rnd_size; i++) {
        rnd[i] = std::atoi(file[index++].c_str());
    }

    InfoLoadResult result;
    result.dimension = dimension;
    result.rnd = rnd;
    return result;
}
