#pragma once

#include <vector>

namespace MLP {

    class Data {
    private:
        int dimension;
        std::vector<double> matrix;
        std::vector<int> rnd;
        std::vector<double> c;
        int rnd_index;

        int to_1D(int i, int j) const {
            return dimension * i + j;
        }

        int to_1D_matrix(int i, int j) const {
            return dimension * i + j;
        }

    public:
        Data();

        void loadData();

        double getDistance(int i, int j) const {
            return matrix[to_1D_matrix(i, j)];
        }

        const std::vector<int> & GetRnd() const { return rnd; }

        int GetDimen() const { return dimension; }

        double GetCost(int i, int j) const {
            return c[to_1D(i, j)];
        }

        int GetRndCrnt() { return rnd[rnd_index++]; }
    };

}
