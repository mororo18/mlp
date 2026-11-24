#include "Data.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <limits>
#include <string>
#include <vector>

#include <sys/resource.h>
#include <sys/time.h>

static double _cputime() {
    struct rusage buf;
    getrusage(RUSAGE_SELF, &buf);
    struct timeval utime = buf.ru_utime;
    struct timeval stime = buf.ru_stime;
    return (utime.tv_sec + utime.tv_usec / 1e6) + (stime.tv_sec + stime.tv_usec / 1e6);
}

struct tData {
    std::vector<double> cost;
    int dimen;

    int T;
    int W;
    int C;

    int REINSERTION;
    int OR_OPT2;
    int OR_OPT3;
    int TWO_OPT;
    int SWAP;

    double EPSILON;

    long reinsert_count;

    std::vector<int> rnd;
    int rnd_index;

    inline double c(int i, int j) const { return cost[i * (dimen + 1) + j]; }
};

struct tSolution {
    std::vector<int> s;   // 1-based access through S(i)
    std::vector<double> seq;
    double cost;
    int stride;

    inline int &S(int i) { return s[i - 1]; }
    inline int S(int i) const { return s[i - 1]; }

    inline double &Q(int t, int j, int i) { return seq[(i * stride + j) * 3 + t]; }
    inline double Q(int t, int j, int i) const { return seq[(i * stride + j) * 3 + t]; }
};

static tSolution solution_init(const tData &data, double cost) {
    tSolution solut;
    solut.s.assign(data.dimen + 1, 0);
    solut.stride = data.dimen + 2;
    solut.seq.assign(3 * solut.stride * solut.stride, 0.0);
    solut.cost = cost;
    return solut;
}

static void update_subseq_info_matrix(tSolution &solut, const tData &data) {
    for (int i = 1; i <= data.dimen + 1; i++) {
        int k = 1 - i - (i == 1);

        solut.Q(data.T, i, i) = 0.0;
        solut.Q(data.C, i, i) = 0.0;
        solut.Q(data.W, i, i) = (i != 1);

        for (int j = i + 1; j <= data.dimen + 1; j++) {
            int j_prev = j - 1;

            solut.Q(data.T, j, i) = data.c(solut.S(j_prev), solut.S(j)) + solut.Q(data.T, j_prev, i);
            solut.Q(data.C, j, i) = solut.Q(data.T, j, i) + solut.Q(data.C, j_prev, i);
            solut.Q(data.W, j, i) = j + k;
        }
    }
    solut.cost = solut.Q(data.C, data.dimen + 1, 1);
}

static int partition(std::vector<int> &arr, int left, int right, const tData &data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j <= right - 1; j++) {
        if (data.c(r, arr[j]) < data.c(r, pivot)) {
            i += 1;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

static void quicksort(std::vector<int> &arr, int left, int right, const tData &data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

static void sort(std::vector<int> &arr, int r, const tData &data) {
    quicksort(arr, 0, (int)arr.size() - 1, data, r);
}

static std::vector<int> construction(double alpha, tData &data) {
    (void)alpha;

    std::vector<int> s;
    s.push_back(1);

    std::vector<int> cList;
    for (int i = 2; i <= data.dimen; i++)
        cList.push_back(i);

    int r = 1;
    while (cList.size() > 0) {
        sort(cList, r, data);

        int r_value = data.rnd[data.rnd_index] + 1;
        data.rnd_index += 1;

        int cN = cList[r_value - 1];

        s.push_back(cN);
        r = cN;
        cList.erase(std::find(cList.begin(), cList.end(), cN));
    }

    s.push_back(1);
    return s;
}

static void swap_nodes(std::vector<int> &s, int i, int j) {
    std::swap(s[i - 1], s[j - 1]);
}

static void reinsert(std::vector<int> &s, int i, int j, int pos) {
    int sz = j - i + 1;

    std::vector<int> blk(s.begin() + (i - 1), s.begin() + j);
    s.insert(s.begin() + (pos - 1), blk.begin(), blk.end());

    if (i < pos) {
        s.erase(s.begin() + (i - 1), s.begin() + (i - 1) + sz);
    } else {
        s.erase(s.begin() + (i + sz - 1), s.begin() + (i + sz - 1) + sz);
    }
}

static bool search_swap(tSolution &solut, const tData &data) {

    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;

    double cost_concat_1 = 0.0;
    double cost_concat_2 = 0.0;
    double cost_concat_3 = 0.0;
    double cost_concat_4 = 0.0;
    double cost_new = 0.0;

    for (int i = 2; i <= data.dimen - 1; i++) {
        int i_prev = i - 1;
        int i_next = i + 1;

        cost_concat_1 =                 solut.Q(data.T, i_prev, 1) + data.c(solut.S(i_prev), solut.S(i_next));
        cost_concat_2 = cost_concat_1 + solut.Q(data.T, i_next, i) + data.c(solut.S(i), solut.S(i_next + 1));

        cost_new = solut.Q(data.C, i_prev, 1)                                                                     +   //       1st subseq
            solut.Q(data.W, i_next, i)                 * (cost_concat_1) + data.c(solut.S(i_next), solut.S(i))     +   // concat 2nd subseq
            solut.Q(data.W, data.dimen + 1, i_next + 1) * (cost_concat_2) + solut.Q(data.C, data.dimen + 1, i_next + 1); // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - data.EPSILON;
            I = i;
            J = i_next;
        }

        for (int j = i_next + 1; j <= data.dimen; j++) {
            int j_prev = j - 1;
            int j_next = j + 1;

            cost_concat_1 =                 solut.Q(data.T, i_prev, 1)      + data.c(solut.S(i_prev), solut.S(j));
            cost_concat_2 = cost_concat_1                                   + data.c(solut.S(j), solut.S(i_next));
            cost_concat_3 = cost_concat_2 + solut.Q(data.T, j_prev, i_next) + data.c(solut.S(j_prev), solut.S(i));
            cost_concat_4 = cost_concat_3                                   + data.c(solut.S(i), solut.S(j_next));

            cost_new = solut.Q(data.C, i_prev, 1)                                                              +   // 1st subseq
                    cost_concat_1 +                                                                                // concat 2nd subseq (single node)
                    solut.Q(data.W, j_prev, i_next)       * cost_concat_2 + solut.Q(data.C, j_prev, i_next)     +   // concat 3rd subseq
                    cost_concat_3 +                                                                                // concat 4th subseq (single node)
                    solut.Q(data.W, data.dimen + 1, j_next) * cost_concat_4 + solut.Q(data.C, data.dimen + 1, j_next); // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - data.EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - data.EPSILON) {
        swap_nodes(solut.s, I, J);
        update_subseq_info_matrix(solut, data);
        return true;
    }

    return false;
}

static bool search_two_opt(tSolution &solut, const tData &data) {
    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;

    double cost_concat_1 = 0.0;
    double cost_concat_2 = 0.0;
    double cost_new = 0.0;

    for (int i = 2; i <= data.dimen - 1; i++) {
        int i_prev = i - 1;
        double rev_seq_cost = solut.Q(data.T, i + 1, i);
        for (int j = i + 2; j <= data.dimen; j++) {
            int j_next = j + 1;

            rev_seq_cost += data.c(solut.S(j - 1), solut.S(j)) * (solut.Q(data.W, j, i) - 1.0);

            cost_concat_1 =                 solut.Q(data.T, i_prev, 1) + data.c(solut.S(j), solut.S(i_prev));
            cost_concat_2 = cost_concat_1 + solut.Q(data.T, j, i)      + data.c(solut.S(j_next), solut.S(i));

            cost_new = solut.Q(data.C, i_prev, 1)                                                                +   //   1st subseq
                    solut.Q(data.W, j, i)                 * cost_concat_1 + rev_seq_cost                         +   // concat 2nd subseq (reversed seq)
                    solut.Q(data.W, data.dimen + 1, j_next) * cost_concat_2 + solut.Q(data.C, data.dimen + 1, j_next); // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - data.EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - data.EPSILON) {
        std::reverse(solut.s.begin() + (I - 1), solut.s.begin() + J);
        update_subseq_info_matrix(solut, data);
        return true;
    }

    return false;
}

static bool search_reinsertion(tSolution &solut, tData &data, int opt) {
    data.reinsert_count += 1;
    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;
    int POS = -1;

    double cost_concat_1 = 0.0;
    double cost_concat_2 = 0.0;
    double cost_concat_3 = 0.0;
    double cost_new = 0.0;

    for (int i = 2; i <= data.dimen - opt + 1; i++) {
        int j = opt + i - 1;
        int i_prev = i - 1;
        int j_next = j + 1;

        for (int k = 1; k <= i_prev - 1; k++) {
            int k_next = k + 1;

            cost_concat_1 =                 solut.Q(data.T, k, 1)             + data.c(solut.S(k), solut.S(i));
            cost_concat_2 = cost_concat_1 + solut.Q(data.T, j, i)             + data.c(solut.S(j), solut.S(k_next));
            cost_concat_3 = cost_concat_2 + solut.Q(data.T, i_prev, k_next)   + data.c(solut.S(i_prev), solut.S(j_next));

            cost_new = solut.Q(data.C, k, 1)                                                                     +   //       1st subseq
                solut.Q(data.W, j, i)                 * cost_concat_1 + solut.Q(data.C, j, i)                    +   // concat 2nd subseq (reinserted seq)
                solut.Q(data.W, i_prev, k_next)       * cost_concat_2 + solut.Q(data.C, i_prev, k_next)          +   // concat 3rd subseq
                solut.Q(data.W, data.dimen + 1, j_next) * cost_concat_3 + solut.Q(data.C, data.dimen + 1, j_next);   // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - data.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (int k = i + opt; k <= data.dimen; k++) {
            int k_next = k + 1;

            cost_concat_1 =                 solut.Q(data.T, i_prev, 1) + data.c(solut.S(i_prev), solut.S(j_next));
            cost_concat_2 = cost_concat_1 + solut.Q(data.T, k, j_next) + data.c(solut.S(k), solut.S(i));
            cost_concat_3 = cost_concat_2 + solut.Q(data.T, j, i)      + data.c(solut.S(j), solut.S(k_next));

            cost_new = solut.Q(data.C, i_prev, 1)                                                                +   //       1st subseq
                solut.Q(data.W, k, j_next)            * cost_concat_1 + solut.Q(data.C, k, j_next)               +   // concat 2nd subseq
                solut.Q(data.W, j, i)                 * cost_concat_2 + solut.Q(data.C, j, i)                    +   // concat 3rd subseq (reinserted seq)
                solut.Q(data.W, data.dimen + 1, k_next) * cost_concat_3 + solut.Q(data.C, data.dimen + 1, k_next);   // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - data.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost - data.EPSILON) {
        reinsert(solut.s, I, J, POS + 1);
        update_subseq_info_matrix(solut, data);
        return true;
    }

    return false;
}

static void RVND(tSolution &solut, tData &data) {
    std::vector<int> neighbd_list = {data.SWAP, data.TWO_OPT, data.REINSERTION, data.OR_OPT2, data.OR_OPT3};

    while (neighbd_list.size() > 0) {

        int i = data.rnd[data.rnd_index] + 1;
        data.rnd_index += 1;

        int neighbd = neighbd_list[i - 1];

        bool improve = false;

        if (neighbd == data.REINSERTION) {
            improve = search_reinsertion(solut, data, data.REINSERTION);
        } else if (neighbd == data.OR_OPT2) {
            improve = search_reinsertion(solut, data, data.OR_OPT2);
        } else if (neighbd == data.OR_OPT3) {
            improve = search_reinsertion(solut, data, data.OR_OPT3);
        } else if (neighbd == data.SWAP) {
            improve = search_swap(solut, data);
        } else if (neighbd == data.TWO_OPT) {
            improve = search_two_opt(solut, data);
        }

        if (improve) {
            neighbd_list = {data.SWAP, data.TWO_OPT, data.REINSERTION, data.OR_OPT2, data.OR_OPT3};
        } else {
            neighbd_list.erase(neighbd_list.begin() + (i - 1));
        }
    }
}

static std::vector<int> perturb(const tSolution &solut_partial, tData &data) {
    std::vector<int> s = solut_partial.s;

    int A_start = 1, A_end = 1;
    int B_start = 1, B_end = 1;

    int size_max = (int)std::floor(s.size() / 10.0);
    size_max = (size_max >= 2 ? size_max : 2);
    int size_min = 2;
    (void)size_max;
    (void)size_min;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = data.rnd[data.rnd_index] + 1;
        data.rnd_index += 1;
        A_end = A_start + data.rnd[data.rnd_index];
        data.rnd_index += 1;

        B_start = data.rnd[data.rnd_index] + 1;
        data.rnd_index += 1;
        B_end = B_start + data.rnd[data.rnd_index];
        data.rnd_index += 1;
    }

    if (A_start < B_start) {
        reinsert(s, B_start, B_end - 1, A_end);
        reinsert(s, A_start, A_end - 1, B_end);
    } else {
        reinsert(s, A_start, A_end - 1, B_end);
        reinsert(s, B_start, B_end - 1, A_end);
    }

    return s;
}

static void GILS_RVND(int Imax, int Iils, const std::vector<double> &R, tData &data, bool verbose) {

    tSolution solut_best = solution_init(data, std::numeric_limits<double>::infinity());
    tSolution solut_partial = solution_init(data, 0.0);
    tSolution solut_crnt = solution_init(data, 0.0);

    for (int i = 1; i <= Imax; i++) {
        int r_value = data.rnd[data.rnd_index] + 1;
        data.rnd_index += 1;

        double alpha = R[r_value - 1];

        if (verbose) {
            std::printf("[+] Local Search %d\n", i);
        }
        solut_crnt.s = construction(alpha, data);
        update_subseq_info_matrix(solut_crnt, data);
        if (verbose) {
            std::printf("\t[+] Constructing Inital Solution.. %.2lf\n", solut_crnt.cost);
        }

        solut_partial.cost = solut_crnt.cost;
        solut_partial.s = solut_crnt.s;

        if (verbose) {
            std::printf("\t[+] Looking for the best Neighbor..\n");
        }
        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(solut_crnt, data);
            if (solut_crnt.cost < solut_partial.cost - data.EPSILON) {
                solut_partial.cost = solut_crnt.cost;
                solut_partial.s = solut_crnt.s;
                iterILS = 0;
            }

            solut_crnt.s = perturb(solut_partial, data);
            update_subseq_info_matrix(solut_crnt, data);

            iterILS += 1;
        }

        update_subseq_info_matrix(solut_partial, data);

        if (solut_partial.cost < solut_best.cost) {
            solut_best.cost = solut_partial.cost;
            solut_best.s = solut_partial.s;
        }

        if (verbose) {
            std::printf("\tCurrent best solution cost: %.2lf\n", solut_best.cost);
        }
    }
    std::printf("COST: %.2lf\n", solut_best.cost);
    std::printf("SOLUTION: ");
    std::printf("[");
    for (size_t i = 0; i < solut_best.s.size(); i++) {
        std::printf("%d", solut_best.s[i]);
        if (i + 1 < solut_best.s.size())
            std::printf(", ");
    }
    std::printf("]\n");
}

int main(int argc, char **argv) {

    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0)
            verbose = true;
    }

    int dimension;
    std::vector<double> cost;
    std::vector<int> rand_values;
    get_instance_info(dimension, cost, rand_values);

    std::vector<double> R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                             0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};

    int Imax = 10;
    int Iils = std::min(dimension, 100);

    tData data;
    data.cost = cost;
    data.dimen = dimension;
    data.T = 0;
    data.W = 1;
    data.C = 2;
    data.REINSERTION = 1;
    data.OR_OPT2 = 2;
    data.OR_OPT3 = 3;
    data.TWO_OPT = 4;
    data.SWAP = 5;
    data.EPSILON = DBL_EPSILON;
    data.reinsert_count = 0;
    data.rnd = rand_values;
    data.rnd_index = 0;

    double _cpu0 = _cputime();
    auto t0 = std::chrono::steady_clock::now();
    GILS_RVND(Imax, Iils, R, data, verbose);
    auto t1 = std::chrono::steady_clock::now();
    double cpu_time = _cputime() - _cpu0;
    double time = std::chrono::duration<double>(t1 - t0).count();

    std::printf("TIME: %.6lf\n", cpu_time);
    std::printf("wall clock (s): %.6lf\n", time);
    std::printf("reinsertion calls %ld\n", data.reinsert_count);

    return 0;
}
