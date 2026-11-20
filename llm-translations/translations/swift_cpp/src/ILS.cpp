#include "ILS.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include "SwiftPrint.hpp"

namespace {

enum class MoveKind {
	swap,
	twoOpt,
	orOpt,
};

struct Move {
	MoveKind kind;
	int blockSize; // only meaningful for orOpt
};

const std::vector<double> R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                              0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                              0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
const std::vector<Move> moves = {{MoveKind::swap, 0},
                                 {MoveKind::twoOpt, 0},
                                 {MoveKind::orOpt, 1},
                                 {MoveKind::orOpt, 2},
                                 {MoveKind::orOpt, 3}};
const double EPSILON = 1e-16;

const double INFINITY_D = std::numeric_limits<double>::infinity();

void reverse(std::vector<int> &self, int i, int j) {
	while (i < j) {
		std::swap(self[i], self[j]);
		i += 1;
		j -= 1;
	}
}

void reinsert(std::vector<int> &self, int blockStart, int blockSize, int at) {
	if (at > blockStart) {
		reverse(self, blockStart, at);
		reverse(self, blockStart, at - blockSize);
		reverse(self, at - blockSize + 1, at);
	} else {
		reverse(self, at + 1, blockStart + blockSize);
		reverse(self, at + 1, at + 1 + blockSize);
		reverse(self, at + 1 + blockSize, blockStart + blockSize);
	}
}

int partition(std::vector<int> &arr, int left, int right, const Data &data, int r) {
	int pivot = arr[right];
	int i = left - 1;

	for (int j = left; j < right; j++) {
		if (data.getDistance(r, arr[j]) < data.getDistance(r, pivot)) {
			i += 1;
			std::swap(arr[i], arr[j]);
		}
	}

	std::swap(arr[i + 1], arr[right]);
	return i + 1;
}

void quickSort(std::vector<int> &arr, int left, int right, const Data &data, int r) {

	if (left >= right) {
		return;
	}

	int pivot = partition(arr, left, right, data, r);
	quickSort(arr, left, pivot - 1, data, r);
	quickSort(arr, pivot + 1, right, data, r);
}

void sort(std::vector<int> &arr, int r, const Data &data) {
	quickSort(arr, 0, (int)arr.size() - 1, data, r);
}

Solution construction(Data &data) {

	Solution s;
	s.sequence = {0};

	(void)data.rndCrnt(); // Mock alpha rnd choice
	std::vector<int> cl;
	for (int i = 1; i < data.dimension; i++) {
		cl.push_back(i);
	}
	int r = 0;
	while (!cl.empty()) {
		sort(cl, r, data);
		int index = data.rndCrnt();

		int c = cl[index];
		r = c;
		s.sequence.push_back(c);
		cl.erase(cl.begin() + index);
	}
	s.sequence.push_back(0);
	s.updateSubseqMatrix(data);

	return s;
}

Solution perturb(Solution s, Data &data) {
	int aStart = 1;
	int aEnd = 1;
	int bStart = 1;
	int bEnd = 1;

	while ((aStart <= bStart && bStart <= aEnd) || (bStart <= aStart && aStart <= bEnd)) {
		aStart = data.rndCrnt();
		aEnd = aStart + data.rndCrnt();

		bStart = data.rndCrnt();
		bEnd = bStart + data.rndCrnt();
	}

	if (bStart < aStart) {
		std::swap(aStart, bStart);
		std::swap(aEnd, bEnd);
	}

	int aBlockSize = aEnd - aStart;
	int bBlockSize = bEnd - bStart;

	reinsert(s.sequence, bStart, bBlockSize, aEnd - 1);
	reinsert(s.sequence, aStart, aBlockSize, bEnd - 1);

	s.updateSubseqMatrix(data);

	return s;
}

bool searchSwap(Solution &s, const Data &data) {

	double costBest = INFINITY_D;
	int bestI = -1;
	int bestJ = -1;

	for (int i = 1; i < data.dimension - 1; i++) {
		int iPrev = i - 1;
		int iNext = i + 1;
		int sIPrev = s.sequence[iPrev];
		int sI = s.sequence[i];
		int sINext = s.sequence[iNext];

		double costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sIPrev, sINext);
		double costConcat2 = costConcat1 + s.getSeq(i, iNext).t + data.getDistance(sI, s.sequence[i + 2]);

		double costNew = s.getSeq(0, iPrev).c
		               + s.getSeq(i, iNext).w                  * costConcat1 + data.getDistance(sINext, sI)
		               + s.getSeq(iNext + 1, data.dimension).w * costConcat2 + s.getSeq(iNext + 1, data.dimension).c;

		if (costNew < costBest) {
			costBest = costNew - EPSILON;
			bestI = i;
			bestJ = iNext;
		}

		for (int j = iNext + 1; j < data.dimension; j++) {
			int jNext = j + 1;
			int jPrev = j - 1;
			int sJPrev = s.sequence[jPrev];
			int sJ = s.sequence[j];
			int sJNext = s.sequence[jNext];

			double costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sIPrev, sJ);
			double costConcat2 = costConcat1 + data.getDistance(sJ, sINext);
			double costConcat3 = costConcat2 + s.getSeq(iNext, jPrev).t + data.getDistance(sJPrev, sI);
			double costConcat4 = costConcat3 + data.getDistance(sI, sJNext);

			double costNew = s.getSeq(0, iPrev).c
			               + costConcat1
			               + s.getSeq(iNext, jPrev).w      * costConcat2 + s.getSeq(iNext, jPrev).c
			               + costConcat3
			               + s.getSeq(jNext, data.dimension).w * costConcat4 + s.getSeq(jNext, data.dimension).c;

			if (costNew < costBest) {
				costBest = costNew - EPSILON;
				bestI = i;
				bestJ = j;
			}
		}
	}

	if (costBest < s.cost) {
		std::swap(s.sequence[bestI], s.sequence[bestJ]);
		s.updateSubseqMatrix(data);
		return true;
	}

	return false;
}

bool searchTwoOpt(Solution &s, const Data &data) {

	double costBest = INFINITY_D;
	int bestI = -1;
	int bestJ = -1;

	for (int i = 1; i < data.dimension - 1; i++) {
		int iPrev = i - 1;
		int sIPrev = s.sequence[iPrev];
		int sI = s.sequence[i];
		double revSeqCost = s.getSeq(i, i + 1).t;

		for (int j = i + 2; j < data.dimension; j++) {
			int jNext = j + 1;
			int sJPrev = s.sequence[j - 1];
			int sJ = s.sequence[j];
			int sJNext = s.sequence[jNext];

			revSeqCost += data.getDistance(sJPrev, sJ) * (s.getSeq(i, j).w - 1.0);

			double costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sJ, sIPrev);
			double costConcat2 = costConcat1 + s.getSeq(i, j).t + data.getDistance(sJNext, sI);

			double costNew = s.getSeq(0, iPrev).c
			               + s.getSeq(i, j).w * costConcat1 + revSeqCost
			               + s.getSeq(jNext, data.dimension).w * costConcat2 + s.getSeq(jNext, data.dimension).c;

			if (costNew < costBest) {
				costBest = costNew - EPSILON;
				bestI = i;
				bestJ = j;
			}
		}
	}

	if (costBest < s.cost) {
		reverse(s.sequence, bestI, bestJ);
		s.updateSubseqMatrix(data);
		return true;
	}

	return false;
}

bool searchOrOpt(Solution &s, int blockSize, const Data &data) {

	double costBest = INFINITY_D;
	int bestI = -1;
	int bestPos = -1;

	for (int i = 1; i <= data.dimension - blockSize; i++) {
		int j = blockSize + i - 1;
		int iPrev = i - 1;
		int jNext = j + 1;
		int sI = s.sequence[i];
		int sJ = s.sequence[j];
		int sIPrev = s.sequence[iPrev];
		int sJNext = s.sequence[jNext];

		for (int k = 0; k < iPrev; k++) {
			int kNext = k + 1;
			int sK = s.sequence[k];
			int sKNext = s.sequence[kNext];

			double costConcat1 = s.getSeq(0, k).t + data.getDistance(sK, sI);
			double costConcat2 = costConcat1 + s.getSeq(i, j).t + data.getDistance(sJ, sKNext);
			double costConcat3 = costConcat2 + s.getSeq(kNext, iPrev).t + data.getDistance(sIPrev, sJNext);

			double costNew = s.getSeq(0, k).c
			               + s.getSeq(i, j).w                  * costConcat1 + s.getSeq(i, j).c
			               + s.getSeq(kNext, iPrev).w          * costConcat2 + s.getSeq(kNext, iPrev).c
			               + s.getSeq(jNext, data.dimension).w * costConcat3 + s.getSeq(jNext, data.dimension).c;

			if (costNew < costBest) {
				costBest = costNew - EPSILON;
				bestI = i;
				bestPos = k;
			}
		}

		for (int k = i + blockSize; k < data.dimension; k++) {
			int kNext = k + 1;
			int sK = s.sequence[k];
			int sKNext = s.sequence[kNext];

			double costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sIPrev, sJNext);
			double costConcat2 = costConcat1 + s.getSeq(jNext, k).t + data.getDistance(sK, sI);
			double costConcat3 = costConcat2 + s.getSeq(i, j).t + data.getDistance(sJ, sKNext);

			double costNew = s.getSeq(0, iPrev).c
			               + s.getSeq(jNext, k).w              * costConcat1 + s.getSeq(jNext, k).c
			               + s.getSeq(i, j).w                  * costConcat2 + s.getSeq(i, j).c
			               + s.getSeq(kNext, data.dimension).w * costConcat3 + s.getSeq(kNext, data.dimension).c;

			if (costNew < costBest) {
				costBest = costNew - EPSILON;
				bestI = i;
				bestPos = k;
			}
		}
	}

	if (costBest < s.cost) {
		reinsert(s.sequence, bestI, blockSize, bestPos);
		s.updateSubseqMatrix(data);
		return true;
	}

	return false;
}

void rvnd(Solution &s, Data &data) {

	std::vector<Move> nl = moves;

	while (!nl.empty()) {
		int index = data.rndCrnt();
		Move move = nl[index];
		bool improved = false;
		switch (move.kind) {
			case MoveKind::swap:   improved = searchSwap(s, data); break;
			case MoveKind::twoOpt: improved = searchTwoOpt(s, data); break;
			case MoveKind::orOpt:  improved = searchOrOpt(s, move.blockSize, data); break;
		}

		if (improved) {
			nl = moves;
			continue;
		}

		nl.erase(nl.begin() + index);
	}
}

} // namespace

void gils_rvnd(int iMax, int iIls, Data &data) {
	Solution bestOfAll;
	bestOfAll.cost = INFINITY_D;

	for (int i = 0; i < iMax; i++) {
		std::cout << "[+] Local Search " << i << "\n";

		Solution localBest = construction(data);
		Solution current = localBest;

		std::cout << "\t[+] Constructing Initial Solution.. \n";
		std::cout << "\t\tSequence: " << describeSequence(current.sequence) << "\n";
		std::cout << "\t\tCost: " << describeDouble(current.cost) << "\n";

		std::cout << "\t[+] Looking for the best Neighbor..\n";
		int iterIls = 0;
		while (iterIls < iIls) {

			rvnd(current, data);

			if (current.cost < localBest.cost) {
				localBest = current;
				iterIls = 0;
			}

			current = perturb(localBest, data);

			iterIls += 1;
		}

		if (localBest.cost < bestOfAll.cost) {
			bestOfAll = localBest;
		}
		std::cout << "\tCurrent Best Cost " << describeDouble(bestOfAll.cost) << "\n";
	}

	std::cout << "COST: " << describeDouble(bestOfAll.cost) << "\n";
}
