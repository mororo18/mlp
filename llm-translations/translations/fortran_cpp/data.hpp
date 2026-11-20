#ifndef _DATA_HPP
#define _DATA_HPP

#include "preproc.hpp"
#include "types.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace data {

    inline void load_matrix(tMatrix &c, std::vector<int> &rnd) {
        std::ifstream io("../distance_matrix");
        std::string line;
        int dimen;
        int value_;
        int rnd_count;

        if (!io.is_open()) {
            std::cout << " could not open ../distance_matrix" << std::endl;
            std::exit(1);
        }

        std::getline(io, line);
        dimen = std::atoi(line.c_str());

        c.allocate(dimen);

        for (int i = 1; i <= dimen; i++) {
            c(i, i) = typeReal(0.0);
            std::getline(io, line);
            std::istringstream values(line);
            for (int j = i + 1; j <= dimen; j++) {
                values >> value_;
                c(i, j) = static_cast<typeReal>(value_);
                c(j, i) = static_cast<typeReal>(value_);
            }
        }

        std::getline(io, line);
        std::getline(io, line);
        std::getline(io, line);
        rnd_count = std::atoi(line.c_str());

        rnd.resize(static_cast<std::size_t>(rnd_count));

        std::cout << " qntd rnd   " << rnd_count << std::endl;

        for (int i = 1; i <= rnd_count; i++) {
            std::getline(io, line);
            rnd[static_cast<std::size_t>(i - 1)] = std::atoi(line.c_str());
        }

        io.close();
    }

}

#endif
