#pragma once

#include <vector>

#include "Data.hpp"

struct SubseqInfo {
	double t = 0.0;
	double c = 0.0;
	double w = 0.0;
};

struct Solution {
	std::vector<int> sequence;
	double cost = 0.0;

	// Subseq matrix
	std::vector<SubseqInfo> seq;
	int dimension = 0;

	// This will copy the subseq info
	// possible variations would be implementing CoW for SubseqInfo
	// or making it a class, benchmark is needed
	const SubseqInfo &getSeq(int i, int j) const {
		return seq[i * dimension + j];
	}

	void updateSubseqMatrix(const Data &data);
};
