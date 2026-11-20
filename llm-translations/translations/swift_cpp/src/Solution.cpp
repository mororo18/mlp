#include "Solution.hpp"

void Solution::updateSubseqMatrix(const Data &data) {
	dimension = (int)sequence.size();

	if (seq.empty()) {
		seq = std::vector<SubseqInfo>(dimension * dimension, SubseqInfo());
	}

	for (int i = 0; i < dimension; i++) {
		int k = 1 - i - (i != 0 ? 0 : 1);

		seq[i * dimension + i].t = 0.0;
		seq[i * dimension + i].c = 0.0;
		seq[i * dimension + i].w = (i != 0 ? 1.0 : 0.0);

		for (int j = i + 1; j < dimension; j++) {
			int jPrev = j - 1;
			seq[i * dimension + j].t = data.getDistance(sequence[jPrev], sequence[j]) + this->getSeq(i, jPrev).t;
			seq[i * dimension + j].c = this->getSeq(i, j).t + this->getSeq(i, jPrev).c;
			seq[i * dimension + j].w = (double)(j + k);
		}
	}

	cost = this->getSeq(0, dimension - 1).c;
}
