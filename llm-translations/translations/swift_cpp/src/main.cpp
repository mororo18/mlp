#include <chrono>
#include <cstdio>
#include <iostream>

#include "Data.hpp"
#include "ILS.hpp"

int main() {
	Data data;

	int iIls = data.dimension < 100 ? data.dimension : 100;

	auto start = std::chrono::steady_clock::now();
	gils_rvnd(10, iIls, data);
	auto end = std::chrono::steady_clock::now();

	double timeInSeconds = std::chrono::duration<double>(end - start).count();

	char buf[64];
	std::snprintf(buf, sizeof(buf), "%.6f", timeInSeconds);
	std::cout << "TIME: " << buf << "\n";

	return 0;
}
