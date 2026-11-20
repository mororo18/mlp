#include "Data.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

namespace {

std::vector<std::string> componentsSeparatedBy(const std::string &s, char sep) {
	std::vector<std::string> out;
	std::string current;
	for (char c : s) {
		if (c == sep) {
			out.push_back(current);
			current.clear();
		} else {
			current += c;
		}
	}
	out.push_back(current);
	return out;
}

} // namespace

Data::Data() {
	std::ifstream file("../distance_matrix");
	std::stringstream buffer;
	buffer << file.rdbuf();
	std::string contents = buffer.str();

	std::vector<std::string> lines = componentsSeparatedBy(contents, '\n');

	int dimension = std::stoi(lines[0]);
	this->dimension = dimension;
	this->matrix = std::vector<double>(dimension * dimension, 0.0);
	for (int i = 0; i <= dimension - 2; i++) {

		this->matrix[i * dimension + i] = 0;
		int j = i + 1;
		std::vector<std::string> values = componentsSeparatedBy(lines[i + 1], ' ');
		values.pop_back();
		for (const std::string &value : values) {

			double converted = std::stod(value);
			this->matrix[i * dimension + j] = converted;
			this->matrix[j * dimension + i] = converted;
			j += 1;
		}
	}

	int rndCount = std::stoi(lines[dimension + 3]);
	for (int i = 0; i < rndCount; i++) {
		int linesIdx = dimension + 4 + i;
		this->rnd.push_back(std::stoi(lines[linesIdx]));
	}
}
