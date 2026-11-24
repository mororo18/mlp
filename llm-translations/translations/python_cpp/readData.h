#ifndef READDATA_H
#define READDATA_H

#include <string>
#include <vector>

typedef std::vector<std::vector<double>> Matrix;

Matrix matrix_fill(int n);

void get_instance_info(int &dimension, Matrix &cost, std::vector<int> &rnd,
                       const std::string &f_name = "../distance_matrix");

#endif
