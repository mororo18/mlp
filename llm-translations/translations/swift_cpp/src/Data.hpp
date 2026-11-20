#pragma once

#include <vector>

// Check if class has better performance over struct
struct Data {
	int dimension = 0;
	std::vector<double> matrix;
	std::vector<int> rnd;
	int rndIdx = 0;

	double getDistance(int i, int j) const {
		return matrix[i * dimension + j];
	}

	int rndCrnt() {
		int value = rnd[rndIdx];
		rndIdx += 1;
		return value;
	}

	Data();
};
