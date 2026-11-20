#pragma once

#include <vector>

#include "Data.hpp"
#include "tSolution.hpp"

namespace MLP {

    class GILS_RVND {
    public:
        static constexpr double EPSILON = 1e-15;

    private:
        Data data;

        int Iils;
        static const int Imax = 10;
        double R[26] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                        0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
        static const int R_size = 26;

        void update_subseq_info_matrix(tSolution & solut);

        void sort(std::vector<int> & arr, int r);
        void Quicksort(std::vector<int> & arr, int left, int right, int r);
        int Partition(std::vector<int> & arr, int left, int right, int r);

        std::vector<int> construction(double alpha);

        void swap(std::vector<int> & s, int i, int j);
        void reverse(std::vector<int> & s, int i, int j);
        void reinsert(std::vector<int> & s, int i, int j, int pos);

        bool search_swap(tSolution & solut);
        bool search_two_opt(tSolution & solut);
        bool search_reinsertion(tSolution & solut, int opt);

        void RVND(tSolution & solut);

        std::vector<int> perturb(const std::vector<int> & sl);

    public:
        GILS_RVND();

        void solve(bool verbose);
    };

}
