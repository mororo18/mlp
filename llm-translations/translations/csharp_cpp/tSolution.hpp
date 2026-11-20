#pragma once

#include <vector>
#include <utility>

namespace MLP {

    class tSolution {
    private:
        std::vector<int> s;
        std::vector<double> seq;
        double cost;
        int d;

    public:
        tSolution(int dimen, double c) {
            seq.assign((size_t)(dimen + 1) * (dimen + 1) * 3, 0.0);
            d = dimen;
            cost = c;
        }

        void StoreSolut(std::vector<int> sl) { s = std::move(sl); }

        inline double GetSeq(int i, int j, int k) const {
            return seq[i * (d + 1) * 3 + j * 3 + k];
        }

        inline void SetSeq(int i, int j, int k, double v) { seq[i * (d + 1) * 3 + j * 3 + k] = v; }

        double GetCost() const { return cost; }
        void SetCost(double v) { cost = v; }

        std::vector<int> & GetSolut() { return s; }
        const std::vector<int> & GetSolut() const { return s; }
        std::vector<int> GetSolutCpy() const { return std::vector<int>(s); }

        int GetPos(int i) const {
            return s[i];
        }
    };

}
