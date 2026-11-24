#ifndef DATA_H
#define DATA_H

#include <vector>

// Reads "../distance_matrix" and returns the instance dimension, the distance
// matrix (1-based indexing, size (dimension+1)x(dimension+1)) and the list of
// pre-generated random values.
void get_instance_info(int &dimension, std::vector<double> &matrix,
                       std::vector<int> &rand_values);

#endif
