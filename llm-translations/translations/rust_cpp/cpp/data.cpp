#include "data.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace data {

LoadResult load(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        std::fprintf(stderr, "Something went wrong reading the file\n");
        std::exit(1);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    const std::string contents = buffer.str();

    std::size_t line_i = 0;

    std::vector<std::string> lines;
    {
        std::istringstream iss(contents);
        std::string l;
        while (std::getline(iss, l)) {
            if (!l.empty() && l.back() == '\r') {
                l.pop_back();
            }
            lines.push_back(l);
        }
    }

    const std::size_t dimension =
        static_cast<std::size_t>(std::stoull(lines[line_i]));
    line_i += 1;

    CostMatrix c(dimension);

    for (std::size_t l_i = 0; l_i < dimension; l_i++) {
        std::istringstream iter(lines[line_i]);
        line_i += 1;

        std::string n;
        bool has = static_cast<bool>(iter >> n);

        const std::size_t i = l_i;
        std::size_t j = l_i + 1;

        while (j < dimension) {
            c(i, j) = std::stod(n);
            c(j, i) = c(i, j);
            has = static_cast<bool>(iter >> n);
            (void)has;
            j += 1;
        }
    }

    line_i += 2;
    const std::size_t rnd_size =
        static_cast<std::size_t>(std::stoull(lines[line_i]));
    line_i += 1;

    std::vector<std::size_t> rnd;
    for (std::size_t k = 0; k < rnd_size; k++) {
        rnd.push_back(static_cast<std::size_t>(std::stoull(lines[line_i])));
        line_i += 1;
    }

    return LoadResult{dimension, std::move(c), std::move(rnd)};
}

} // namespace data
