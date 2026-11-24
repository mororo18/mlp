#include "readData.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

typedef std::array<double, 3> tSubseq;
typedef std::vector<std::vector<tSubseq>> tSeq;

// reproduces Python's repr() of a float, so the output is identical
static std::string py_float(double v) {
    if (std::isinf(v))
        return v > 0 ? "inf" : "-inf";
    if (std::isnan(v))
        return "nan";

    char buf[64];
    for (int p = 1; p <= 17; p++) {
        std::snprintf(buf, sizeof(buf), "%.*g", p, v);
        if (std::strtod(buf, nullptr) == v)
            break;
    }

    std::string s(buf);
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos)
        s += ".0";
    return s;
}

class tInfo {
public:
    int IT = 0;
    int dimen;
    Matrix cost;

    static const int T = 0;
    static const int C = 1;
    static const int W = 2;

    static const int SWAP = 0;
    static const int REINSERTION = 1;
    static const int OR_OPT_2 = 2;
    static const int OR_OPT_3 = 3;
    static const int TWO_OPT = 4;

    static constexpr double EPSILON = 1e-13;

    std::vector<int> rnd;
    int rnd_index = 0;

    tInfo(int dimension, const Matrix &cost, const std::vector<int> &rnd)
        : dimen(dimension), cost(cost), rnd(rnd) {}
};

class tSolution {
public:
    std::vector<int> s;
    tSeq seq;
    double cost;

    tSolution(const std::vector<int> &s, const tSeq &seq, double cost)
        : s(s), seq(seq), cost(cost) {}
};

tSeq subseq_fill(const tInfo &info) {
    tSeq arr;

    for (int i = 0; i < info.dimen + 1; i++) {
        arr.push_back(std::vector<tSubseq>());
        for (int j = 0; j < info.dimen + 1; j++) {
            arr[i].push_back({0.0, 0.0, 0.0});
        }
    }

    return arr;
}

int partition(std::vector<int> &arr, int left, int right, const tInfo &info, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (info.cost[r][arr[j]] < info.cost[r][pivot]) {
            i += 1;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

void quicksort(std::vector<int> &arr, int left, int right, const tInfo &info, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, info, r);
        quicksort(arr, left, pivot - 1, info, r);
        quicksort(arr, pivot + 1, right, info, r);
    }
}

void sort(std::vector<int> &arr, int r, const tInfo &info) {
    quicksort(arr, 0, (int)arr.size() - 1, info, r);
}

std::vector<int> construction(double alpha, tInfo &info) {
    std::vector<int> s = {0};
    std::vector<int> c_list;
    for (int i = 1; i < info.dimen; i++)
        c_list.push_back(i);

    int r = 0;
    while (c_list.size() > 0) {

        int i = (int)(c_list.size() * alpha) + 1;
        (void)i;

        sort(c_list, r, info);

        int r_value = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        int c = c_list[r_value];

        s.push_back(c);
        r = c;
        c_list.erase(std::find(c_list.begin(), c_list.end(), c));
    }

    s.push_back(0);

    return s;
}

void update_subseq_info_matrix(tSolution &solut, const tInfo &info) {
    const int T = tInfo::T;
    const int C = tInfo::C;
    const int W = tInfo::W;
    const int dimen = info.dimen;
    const Matrix &cost = info.cost;

    for (int i = 0; i < dimen + 1; i++) {
        int k = 1 - i - (int)(!i);

        solut.seq[i][i][T] = 0.0;
        solut.seq[i][i][C] = 0.0;
        solut.seq[i][i][W] = (double)(int)(!(i == 0));

        for (int j = i + 1; j < dimen + 1; j++) {
            int j_prev = j - 1;

            solut.seq[i][j][T] = cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][T];
            solut.seq[i][j][C] = solut.seq[i][j][T] + solut.seq[i][j_prev][C];
            solut.seq[i][j][W] = (double)(j + k);
        }
    }

    solut.cost = solut.seq[0][dimen][C] - tInfo::EPSILON;
}

void swap(std::vector<int> &s, int i, int j) {
    std::swap(s[i], s[j]);
}

void reverse(std::vector<int> &s, int i, int j) {
    std::reverse(s.begin() + i, s.begin() + j + 1);
}

void reinsert(std::vector<int> &s, int i, int j, int pos) {
    std::vector<int> sub(s.begin() + i, s.begin() + j + 1);
    if (i < pos) {
        s.insert(s.begin() + pos, sub.begin(), sub.end());
        s.erase(s.begin() + i, s.begin() + j + 1);
    } else {
        s.erase(s.begin() + i, s.begin() + j + 1);
        s.insert(s.begin() + pos, sub.begin(), sub.end());
    }
}

bool search_swap(tSolution &solut, const tInfo &info) {
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    double cost_concat_4;

    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;

    const int dimen = info.dimen;
    const Matrix &cost = info.cost;
    const std::vector<int> &s = solut.s;
    const tSeq &seq = solut.seq;

    const int T = tInfo::T;
    const int C = tInfo::C;
    const int W = tInfo::W;

    for (int i = 1; i < dimen - 1; i++) {
        int i_prev = i - 1;
        int i_next = i + 1;

        cost_concat_1 = seq[0][i_prev][T] + cost[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[i][i_next][T] + cost[s[i]][s[i_next + 1]];

        cost_new = seq[0][i_prev][C]                                            +
                seq[i][i_next][W]   * (cost_concat_1) + cost[s[i_next]][s[i]]   +
                seq[i_next + 1][dimen][W] * (cost_concat_2) + seq[i_next + 1][dimen][C];

        if (cost_new < cost_best) {
            cost_best = cost_new - tInfo::EPSILON;
            I = i;
            J = i_next;
        }

        for (int j = i_next + 1; j < dimen; j++) {
            int j_next = j + 1;
            int j_prev = j - 1;

            cost_concat_1 = seq[0][i_prev][T] + cost[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1 + cost[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev][T] + cost[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3 + cost[s[i]][s[j_next]];

            cost_new = seq[0][i_prev][C]                                                +
                    cost_concat_1                                                       +
                    seq[i_next][j_prev][W] * cost_concat_2 + seq[i_next][j_prev][C]     +
                    cost_concat_3                                                       +
                    seq[j_next][dimen][W] * cost_concat_4 + seq[j_next][dimen][C];

            if (cost_new < cost_best) {
                cost_best = cost_new - tInfo::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - tInfo::EPSILON) {
        swap(solut.s, I, J);
        update_subseq_info_matrix(solut, info);
        return true;
    }

    return false;
}

bool search_two_opt(tSolution &solut, const tInfo &info) {
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;

    const int dimen = info.dimen;
    const Matrix &cost = info.cost;
    const std::vector<int> &s = solut.s;
    const tSeq &seq = solut.seq;

    const int T = tInfo::T;
    const int C = tInfo::C;
    const int W = tInfo::W;

    for (int i = 1; i < dimen - 1; i++) {
        int i_prev = i - 1;

        double reverse_cost = seq[i][i + 1][T];
        for (int j = i + 2; j < dimen; j++) {
            int j_next = j + 1;

            reverse_cost += cost[s[j - 1]][s[j]] * (seq[i][j][W] - 1);

            cost_concat_1 = seq[0][i_prev][T] + cost[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[i][j][T] + cost[s[j_next]][s[i]];

            cost_new = seq[0][i_prev][C]                                        +
                    seq[i][j][W] * cost_concat_1 + reverse_cost                 +
                    seq[j_next][dimen][W] * (cost_concat_2) + seq[j_next][dimen][C];

            if (cost_new < cost_best) {
                cost_best = cost_new - tInfo::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - tInfo::EPSILON) {
        reverse(solut.s, I, J);
        update_subseq_info_matrix(solut, info);
        return true;
    }

    return false;
}

bool search_reinsertion(tSolution &solut, const tInfo &info, int opt, const tSeq &seq) {
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;

    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;
    int POS = -1;
    const int MAX = info.dimen - opt;

    const int dimen = info.dimen;
    const Matrix &cost = info.cost;
    const std::vector<int> &s = solut.s;
    const int T = tInfo::T;
    const int C = tInfo::C;
    const int W = tInfo::W;

    for (int i = 1; i < MAX + 1; i++) {
        int j = opt + i - 1;

        int j_next = j + 1;
        int i_prev = i - 1;

        for (int k = 0; k < i_prev; k++) {
            int k_next = k + 1;

            cost_concat_1 =                 seq[0]     [k]     [T] + cost[s[k]]     [s[i]];
            cost_concat_2 = cost_concat_1 + seq[i]     [j]     [T] + cost[s[j]]     [s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev][T] + cost[s[i_prev]][s[j_next]];

            cost_new = seq[0][k][C]                                                     +
                    seq[i]     [j]     [W]  * cost_concat_1 + seq[i]     [j]     [C]    +
                    seq[k_next][i_prev][W]  * cost_concat_2 + seq[k_next][i_prev][C]    +
                    seq[j_next][dimen] [W]  * cost_concat_3 + seq[j_next][dimen] [C];

            if (cost_new < cost_best) {
                cost_best = cost_new - tInfo::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (int k = i + opt; k < dimen; k++) {
            int k_next = k + 1;

            cost_concat_1 =                 seq[0]     [i_prev][T] + cost[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[j_next][k]     [T] + cost[s[k]]     [s[i]];
            cost_concat_3 = cost_concat_2 + seq[i]     [j]     [T] + cost[s[j]]     [s[k_next]];

            cost_new = seq[0][i_prev][C]                                            +
                    seq[j_next][k]    [W]   * cost_concat_1 + seq[j_next][k]    [C] +
                    seq[i]     [j]    [W]   * cost_concat_2 + seq[i]     [j]    [C] +
                    seq[k_next][dimen][W]   * cost_concat_3 + seq[k_next][dimen][C];

            if (cost_new < cost_best) {
                cost_best = cost_new - tInfo::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost - tInfo::EPSILON) {
        reinsert(solut.s, I, J, POS + 1);
        update_subseq_info_matrix(solut, info);
        return true;
    }

    return false;
}

void RVND(tSolution &solut, tInfo &info) {

    std::vector<int> neighbd_list = {tInfo::SWAP, tInfo::TWO_OPT, tInfo::REINSERTION,
                                     tInfo::OR_OPT_2, tInfo::OR_OPT_3};

    while (neighbd_list.size() > 0) {
        int i = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        int neighbd = neighbd_list[i];

        bool improve = false;

        if (neighbd == tInfo::SWAP) {
            improve = search_swap(solut, info);
        } else if (neighbd == tInfo::TWO_OPT) {
            improve = search_two_opt(solut, info);
        } else if (neighbd == tInfo::REINSERTION) {
            improve = search_reinsertion(solut, info, tInfo::REINSERTION, solut.seq);
        } else if (neighbd == tInfo::OR_OPT_2) {
            improve = search_reinsertion(solut, info, tInfo::OR_OPT_2, solut.seq);
        } else if (neighbd == tInfo::OR_OPT_3) {
            improve = search_reinsertion(solut, info, tInfo::OR_OPT_3, solut.seq);
        }

        if (improve == true) {
            neighbd_list = {tInfo::SWAP, tInfo::TWO_OPT, tInfo::REINSERTION,
                            tInfo::OR_OPT_2, tInfo::OR_OPT_3};

        } else {
            neighbd_list.erase(neighbd_list.begin() + i);
        }
    }
}

std::vector<int> perturb(const std::vector<int> &sl, tInfo &info) {
    std::vector<int> s = sl;

    int A_start = 1, A_end = 1;
    int B_start = 1, B_end = 1;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        A_end = A_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;

        B_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        B_end = B_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;
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

void GILS_RVND(int Imax, int Iils, const std::vector<double> &R, tInfo &info,
               bool verbose = false) {

    tSolution solut_partial(std::vector<int>(info.dimen + 1, 0), subseq_fill(info), 0.0);
    tSolution solut_crnt(std::vector<int>(info.dimen + 1, 0), subseq_fill(info), 0.0);
    tSolution solut_best(std::vector<int>(info.dimen + 1, 0), subseq_fill(info),
                         std::numeric_limits<double>::infinity());

    for (int i = 0; i < Imax; i++) {
        int r_value = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        double alpha = R[r_value];

        if (verbose)
            std::cout << "[+] Local Search " << i + 1 << std::endl;
        solut_crnt.s = construction(alpha, info);
        update_subseq_info_matrix(solut_crnt, info);

        solut_partial.s = solut_crnt.s;
        solut_partial.cost = solut_crnt.cost;
        if (verbose)
            std::cout << "\t[+] Constructing Inital Solution.. " << py_float(solut_crnt.cost)
                      << std::endl;

        if (verbose)
            std::cout << "\t[+] Looking for the best Neighbor.." << std::endl;
        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(solut_crnt, info);

            if (solut_crnt.cost < solut_partial.cost) {
                solut_partial.s = solut_crnt.s;
                solut_partial.cost = solut_crnt.cost;
                iterILS = 0;
            }

            solut_crnt.s = perturb(solut_partial.s, info);
            update_subseq_info_matrix(solut_crnt, info);
            iterILS += 1;
        }

        update_subseq_info_matrix(solut_partial, info);

        if (solut_partial.cost < solut_best.cost) {
            solut_best.s = solut_partial.s;
            solut_best.cost = solut_partial.cost;
        }

        if (verbose)
            std::cout << "\tCurrent best solution cost: " << py_float(solut_best.cost) << std::endl;
    }

    std::cout << "COST: " << py_float(solut_best.cost) << std::endl;

    std::cout << "SOLUTION: [";
    for (size_t i = 0; i < solut_best.s.size(); i++) {
        if (i > 0)
            std::cout << ", ";
        std::cout << solut_best.s[i];
    }
    std::cout << "]" << std::endl;
}

int main(int argc, char **argv) {
    int dimension;
    Matrix cost;
    std::vector<int> rnd;

    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-v" || std::string(argv[i]) == "--verbose")
            verbose = true;
    }

    get_instance_info(dimension, cost, rnd);

    std::vector<double> R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
                             0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                             0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};

    int Imax = 10;
    int Iils = std::min(dimension, 100);

    tInfo info(dimension, cost, rnd);

    clock_t start = clock();
    GILS_RVND(Imax, Iils, R, info, verbose);
    std::cout << "TIME: " << py_float((double)(clock() - start) / CLOCKS_PER_SEC) << " "
              << std::endl;

    std::cout << "ITERACOES:  " << info.IT << std::endl;

    return 0;
}
