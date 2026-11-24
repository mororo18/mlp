#ifndef DATA_H
#define DATA_H

#include <vector>

struct InfoLoadResult {
    int dimension;
    std::vector<int> rnd;
};

InfoLoadResult info_load(std::vector<std::vector<double> > & c);

#endif
