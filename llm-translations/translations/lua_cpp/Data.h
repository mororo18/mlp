#ifndef DATA_H
#define DATA_H

#include <vector>

// Mirrors Data.lua: reads "../distance_matrix", fills the 1-based cost
// matrix `c` and the 1-based random-number list `rnd`.
// Returns the dimension; `rnd_size` receives the amount of rnd entries.
int readData(std::vector<std::vector<double>> &c, std::vector<int> &rnd, int &rnd_size);

#endif
