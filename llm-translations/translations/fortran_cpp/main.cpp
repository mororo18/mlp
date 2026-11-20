#include "preproc.hpp"

#ifdef DEBUG
#include "assert.hpp"
#endif

#include "data.hpp"
#include "types.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

void print_matrix(const tMatrix &c) {
    int dimen;
    int i;
    int j;

    dimen = c.n;

    for (i = 1; i <= 2; i++) {
        for (j = 1; j <= dimen; j++) {
            std::printf("%8.3f ", static_cast<double>(c(i, j)));
        }
        std::printf("\n\n");
        std::printf("\n");
    }
}

void print_info(const tInfo &info) {
    std::cout << " " << info.dimen << std::endl;
    print_matrix(info.cost);
}

void subseq_load(tSolution &sol, const tInfo &info) {
    int i;
    int j;
    int k;
    int j_prev;

    for (i = 1; i <= info.dimen + 1; i++) {
        k = 1 - i - (i != 1 ? 0 : 1);

        sol.seq_at(i, i).T = typeReal(0.0);
        sol.seq_at(i, i).C = typeReal(0.0);
        sol.seq_at(i, i).W = (i != 1 ? typeReal(1.0) : typeReal(0.0));

        for (j = i + 1; j <= info.dimen + 1; j++) {

            j_prev = j - 1;

            sol.seq_at(j, i).T = info.cost(sol.s_at(j_prev), sol.s_at(j)) + sol.seq_at(j_prev, i).T;
            sol.seq_at(j, i).C = sol.seq_at(j, i).T + sol.seq_at(j_prev, i).C;
            sol.seq_at(j, i).W = static_cast<typeReal>(j + k);

        }
    }

    sol.cost = sol.seq_at(info.dimen + 1, 1).C;
}

int partition_(std::vector<int> &arr, int left, int right, const tInfo &info, int r) {
    int pivot, i, j, temp;

    pivot = arr[right - 1];
    i = left - 1;

    for (j = left; j <= right - 1; j++) {
        if (info.cost(r, arr[j - 1]) < info.cost(r, pivot)) {
            i = i + 1;
            temp = arr[i - 1];
            arr[i - 1] = arr[j - 1];
            arr[j - 1] = temp;
        }
    }
    temp = arr[i + 1 - 1];
    arr[i + 1 - 1] = arr[right - 1];
    arr[right - 1] = temp;

    return i + 1;
}

void quicksort(std::vector<int> &arr, int left, int right, const tInfo &info, int r) {
    int pivot;

    if (left < right) {
        pivot = partition_(arr, left, right, info, r);
        quicksort(arr, left, pivot - 1, info, r);
        quicksort(arr, pivot + 1, right, info, r);
    }
}

void sort(std::vector<int> &arr, int len, int r, const tInfo &info) {
    quicksort(arr, 1, len, info, r);
}

void arr_shift(std::vector<int> &arr, int src, int tgt, int sz) {
    int i;
    int diff;

    if (src < tgt) {
        diff = tgt - src;
        for (i = src + sz - 1; i >= src; i--) {
            arr[i + diff - 1] = arr[i - 1];
        }
    } else {
        for (i = 1; i <= sz; i++) {
            arr[tgt + i - 1 - 1] = arr[i + src - 1 - 1];
        }
    }
}

std::vector<int> construction(typeReal alpha, tInfo &info) {
    (void)alpha;

    int i;
    int j;
    std::vector<int> s(static_cast<std::size_t>(info.dimen + 1));
    std::vector<int> cL(static_cast<std::size_t>(info.dimen - 1));
    int cL_size;
    int cN;
    int r;
    int index_;

    cL_size = info.dimen - 1;

    for (i = 2; i <= info.dimen; i++) {
        cL[static_cast<std::size_t>(i - 2)] = i;
    }

    s[0] = 1;
    r = 1;
    for (j = 2; j <= info.dimen; j++) {
        sort(cL, cL_size, r, info);

        index_ = info.rnd_at(info.rnd_index) + 1;
        info.rnd_index = info.rnd_index + 1;

        cN = cL[static_cast<std::size_t>(index_ - 1)];

        arr_shift(cL, index_ + 1, index_, cL_size - index_);
        cL_size = cL_size - 1;

        s[static_cast<std::size_t>(j - 1)] = cN;
        r = cN;
    }

    s[static_cast<std::size_t>(info.dimen)] = 1;

    return s;
}

void swap(std::vector<int> &s, int i, int j) {
    int tmp;

    tmp = s[i - 1];
    s[i - 1] = s[j - 1];
    s[j - 1] = tmp;
}

void reverse(std::vector<int> &s, int i, int j) {
    int index_;
    int bound;

    bound = j;
    for (index_ = i; index_ <= (i + j) / 2; index_++) {
        swap(s, index_, bound);
        bound = bound - 1;
    }
}

void reinsert(std::vector<int> &s, int i, int j, int pos) {
    std::vector<int> sub(static_cast<std::size_t>(j - i + 1));
    int sz;

    for (int k = i; k <= j; k++) {
        sub[static_cast<std::size_t>(k - i)] = s[k - 1];
    }

    if (i < pos) {
        sz = pos - j - 1;
        arr_shift(s, j + 1, i, sz);
        for (int k = i + sz; k <= pos - 1; k++) {
            s[k - 1] = sub[static_cast<std::size_t>(k - (i + sz))];
        }
    } else {
        sz = i - pos;
        arr_shift(s, pos, j + 1 - sz, sz);
        for (int k = pos; k <= pos + j - i; k++) {
            s[k - 1] = sub[static_cast<std::size_t>(k - pos)];
        }
    }
}

void is_feasible(const tSolution &solut, bool &ret) {
    int i;
    int n;

    n = solut.s_size - 1;
    std::vector<bool> check(static_cast<std::size_t>(n), false);

    for (i = 1; i <= n; i++) {
        check[static_cast<std::size_t>(solut.s_at(i) - 1)] = true;
    }

    ret = true;
    for (i = 1; i <= n; i++) {
        if (check[static_cast<std::size_t>(i - 1)] == false) {
            ret = false;
        }
    }
}

void search_swap(tSolution &solut, const tInfo &info, bool &ret) {
    int i;
    int j;
    int i_prev;
    int i_next;
    int j_prev;
    int j_next;

    int I_best = 0;
    int J_best = 0;

    typeReal cost_best;
    typeReal cost_new;
    typeReal cost_concat_1;
    typeReal cost_concat_2;
    typeReal cost_concat_3;
    typeReal cost_concat_4;

#ifdef DEBUG
    bool fsb;
#endif

    cost_best = info.fmax;
    cost_new = typeReal(0.0);
    cost_concat_1 = typeReal(0.0);
    cost_concat_2 = typeReal(0.0);
    cost_concat_3 = typeReal(0.0);
    cost_concat_4 = typeReal(0.0);

    for (i = 2; i <= info.dimen - 1; i++) {
        i_prev = i - 1;
        i_next = i + 1;

        cost_concat_1 =                 solut.seq_at(i_prev, 1).T + info.cost(solut.s_at(i_prev), solut.s_at(i_next));
        cost_concat_2 = cost_concat_1 + solut.seq_at(i_next, i).T + info.cost(solut.s_at(i), solut.s_at(i_next + 1));

        cost_new = solut.seq_at(i_prev, 1).C                                                                             +
                solut.seq_at(i_next, i).W                * (cost_concat_1) + info.cost(solut.s_at(i_next), solut.s_at(i)) +
                solut.seq_at(info.dimen + 1, i_next + 1).W * (cost_concat_2) + solut.seq_at(info.dimen + 1, i_next + 1).C;

        if (cost_new < cost_best) {
            cost_best = cost_new;
            I_best = i;
            J_best = i_next;
        }

        for (j = i_next + 1; j <= info.dimen; j++) {

            j_prev = j - 1;
            j_next = j + 1;

            cost_concat_1 =                 solut.seq_at(i_prev, 1).T      + info.cost(solut.s_at(i_prev), solut.s_at(j));
            cost_concat_2 = cost_concat_1                                  + info.cost(solut.s_at(j), solut.s_at(i_next));
            cost_concat_3 = cost_concat_2 + solut.seq_at(j_prev, i_next).T + info.cost(solut.s_at(j_prev), solut.s_at(i));
            cost_concat_4 = cost_concat_3                                  + info.cost(solut.s_at(i), solut.s_at(j_next));

            cost_new = solut.seq_at(i_prev, 1).C                                                                    +   // 1st subseq
                    cost_concat_1 +                                                                                     // concat 2nd subseq (single node)
                    solut.seq_at(j_prev, i_next).W       * cost_concat_2 + solut.seq_at(j_prev, i_next).C           +   // concat 3rd subseq
                    cost_concat_3 +                                                                                     // concat 4th subseq (single node)
                    solut.seq_at(info.dimen + 1, j_next).W * cost_concat_4 + solut.seq_at(info.dimen + 1, j_next).C;    // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I_best = i;
                J_best = j;
            }

        }

    }

    if (cost_best < solut.cost - std::numeric_limits<float>::epsilon()) {
        swap(solut.s, I_best, J_best);
        subseq_load(solut, info);
#ifdef DEBUG
        fsb = false;
        is_feasible(solut, fsb);
        ASSERT(fsb, true);
        ASSERT(cost_best, solut.cost);
#endif
        ret = true;
    } else {
        ret = false;
    }
}

void search_two_opt(tSolution &solut, const tInfo &info, bool &ret) {
    int i;
    int j;
    int i_prev;
    int j_next;

    int I_best = 0;
    int J_best = 0;

    typeReal cost_best;
    typeReal cost_new;
    typeReal cost_concat_1;
    typeReal cost_concat_2;
    typeReal rev_seq_cost;

#ifdef DEBUG
    bool fsb;
#endif

    cost_best = info.fmax;
    cost_new = typeReal(0.0);
    cost_concat_1 = typeReal(0.0);
    cost_concat_2 = typeReal(0.0);

    for (i = 2; i <= info.dimen - 1; i++) {
        i_prev = i - 1;
        rev_seq_cost = solut.seq_at(i + 1, i).T;

        for (j = i + 2; j <= info.dimen; j++) {
            j_next = j + 1;

            rev_seq_cost = rev_seq_cost + info.cost(solut.s_at(j - 1), solut.s_at(j)) * (solut.seq_at(j, i).W - typeReal(1.0));

            cost_concat_1 =                 solut.seq_at(i_prev, 1).T + info.cost(solut.s_at(j), solut.s_at(i_prev));
            cost_concat_2 = cost_concat_1 + solut.seq_at(j, i).T      + info.cost(solut.s_at(j_next), solut.s_at(i));

            cost_new = solut.seq_at(i_prev, 1).C                                                                  +   //   1st subseq
                    solut.seq_at(j, i).W                   * cost_concat_1 + rev_seq_cost                         +   // concat 2nd subseq (reversed seq)
                    solut.seq_at(info.dimen + 1, j_next).W * cost_concat_2 + solut.seq_at(info.dimen + 1, j_next).C;  // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I_best = i;
                J_best = j;
            }

        }

    }

    if (cost_best < solut.cost) {
        reverse(solut.s, I_best, J_best);
        subseq_load(solut, info);
#ifdef DEBUG
        fsb = false;
        is_feasible(solut, fsb);
        ASSERT(fsb, true);
        ASSERT(cost_best, solut.cost);
#endif
        ret = true;
    } else {
        ret = false;
    }
}

void search_reinsertion(tSolution &solut, tInfo &info, int opt, bool &ret) {
    int i;
    int j;
    int k;
    int i_prev;
    int j_next;
    int k_next;

    int I_best = 0;
    int J_best = 0;
    int POS_best = 0;

    typeReal cost_best;
    typeReal cost_new;
    typeReal cost_concat_1;
    typeReal cost_concat_2;
    typeReal cost_concat_3;

#ifdef DEBUG
    bool fsb;
#endif

    info.reinsert_call = info.reinsert_call + 1;

    cost_best = info.fmax;
    cost_new = typeReal(0.0);
    cost_concat_1 = typeReal(0.0);
    cost_concat_2 = typeReal(0.0);
    cost_concat_3 = typeReal(0.0);

    for (i = 2; i <= info.dimen - opt + 1; i++) {
        j = opt + i - 1;
        i_prev = i - 1;
        j_next = j + 1;

        for (k = 1; k <= i_prev - 1; k++) {
            k_next = k + 1;

            cost_concat_1 =                 solut.seq_at(k, 1).T           + info.cost(solut.s_at(k), solut.s_at(i));
            cost_concat_2 = cost_concat_1 + solut.seq_at(j, i).T           + info.cost(solut.s_at(j), solut.s_at(k_next));
            cost_concat_3 = cost_concat_2 + solut.seq_at(i_prev, k_next).T + info.cost(solut.s_at(i_prev), solut.s_at(j_next));

            cost_new = solut.seq_at(k, 1).C                                                                          +   //       1st subseq
                    solut.seq_at(j, i).W                   * cost_concat_1 + solut.seq_at(j, i).C                    +   // concat 2nd subseq (reinserted seq)
                    solut.seq_at(i_prev, k_next).W         * cost_concat_2 + solut.seq_at(i_prev, k_next).C          +   // concat 3rd subseq
                    solut.seq_at(info.dimen + 1, j_next).W * cost_concat_3 + solut.seq_at(info.dimen + 1, j_next).C;     // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I_best = i;
                J_best = j;
                POS_best = k;
            }

        }

        for (k = i + opt; k <= info.dimen; k++) {
            k_next = k + 1;

            cost_concat_1 =                 solut.seq_at(i_prev, 1).T  + info.cost(solut.s_at(i_prev), solut.s_at(j_next));
            cost_concat_2 = cost_concat_1 + solut.seq_at(k, j_next).T  + info.cost(solut.s_at(k), solut.s_at(i));
            cost_concat_3 = cost_concat_2 + solut.seq_at(j, i).T       + info.cost(solut.s_at(j), solut.s_at(k_next));

            cost_new = solut.seq_at(i_prev, 1).C                                                                     +   //       1st subseq
                    solut.seq_at(k, j_next).W              * cost_concat_1 + solut.seq_at(k, j_next).C               +   // concat 2nd subseq
                    solut.seq_at(j, i).W                   * cost_concat_2 + solut.seq_at(j, i).C                    +   // concat 3rd subseq (reinserted seq)
                    solut.seq_at(info.dimen + 1, k_next).W * cost_concat_3 + solut.seq_at(info.dimen + 1, k_next).C;     // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I_best = i;
                J_best = j;
                POS_best = k;
            }
        }
    }

    if (cost_best < solut.cost) {
        reinsert(solut.s, I_best, J_best, POS_best + 1);
        subseq_load(solut, info);

#ifdef DEBUG
        fsb = false;
        is_feasible(solut, fsb);
        ASSERT(fsb, true);
        ASSERT(cost_best, solut.cost);
#endif
        ret = true;
    } else {
        ret = false;
    }
}

void RVND(tSolution &sol, tInfo &info, int &it) {
    std::vector<int> neighbd_list(5);
    int nl_size;
    int neighbd;

    int index_;
    int total;
    bool improve_flag;
#ifdef DEBUG
    bool fsb;
#endif
    bool yes;

    std::vector<int> improv(5, 0);

    nl_size = 5;
    neighbd_list[0] = info.SWAP;
    neighbd_list[1] = info.TWO_OPT;
    neighbd_list[2] = info.REINSERTION;
    neighbd_list[3] = info.OR_OPT_2;
    neighbd_list[4] = info.OR_OPT_3;
    yes = false;
    total = 0;
    while (nl_size > 0) {
        index_ = info.rnd_at(info.rnd_index) + 1;
        info.rnd_index = info.rnd_index + 1;

#ifdef DEBUG
        ASSERT(index_ > nl_size, false);
#endif

        neighbd = neighbd_list[static_cast<std::size_t>(index_ - 1)];

        improve_flag = false;
        if (neighbd == info.REINSERTION) {
            search_reinsertion(sol, info, info.REINSERTION, improve_flag);
        } else if (neighbd == info.OR_OPT_2) {
            search_reinsertion(sol, info, info.OR_OPT_2, improve_flag);
        } else if (neighbd == info.OR_OPT_3) {
            search_reinsertion(sol, info, info.OR_OPT_3, improve_flag);
        } else if (neighbd == info.SWAP) {
            search_swap(sol, info, improve_flag);
        } else if (neighbd == info.TWO_OPT) {
            search_two_opt(sol, info, improve_flag);
        }

#ifdef DEBUG
        // feasibility check
        fsb = false;
        is_feasible(sol, fsb);
        ASSERT(fsb, true);
#endif

        if (improve_flag == true) {
            neighbd_list[0] = info.SWAP;
            neighbd_list[1] = info.TWO_OPT;
            neighbd_list[2] = info.REINSERTION;
            neighbd_list[3] = info.OR_OPT_2;
            neighbd_list[4] = info.OR_OPT_3;
            nl_size = 5;
            improv[static_cast<std::size_t>(neighbd - 1)] = improv[static_cast<std::size_t>(neighbd - 1)] + 1;
            yes = true;

        } else {
            arr_shift(neighbd_list, index_ + 1, index_, nl_size - index_);
            nl_size = nl_size - 1;
            total = total + 1;
        }

        it = it + 1;

    }

    (void)yes;
    (void)total;
}

void perturb(tSolution &solut_crnt, const tSolution &solut_part, tInfo &info) {
    int A_start, A_end;
    int B_start, B_end;
    int size_max;
    int size_min;

    A_start = 1;
    A_end = 1;
    B_start = 1;
    B_end = 1;

    solut_crnt.s = solut_part.s;

    size_max = (info.dimen + 1) / 10;
    size_max = (size_max >= 2 ? size_max : 2);
    size_min = 2;

    (void)size_max;
    (void)size_min;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        A_start = info.rnd_at(info.rnd_index) + 1;
        info.rnd_index = info.rnd_index + 1;
        A_end = A_start + info.rnd_at(info.rnd_index);
        info.rnd_index = info.rnd_index + 1;

        B_start = info.rnd_at(info.rnd_index) + 1;
        info.rnd_index = info.rnd_index + 1;
        B_end = B_start + info.rnd_at(info.rnd_index);
        info.rnd_index = info.rnd_index + 1;
    }

    if (A_start < B_start) {
        reinsert(solut_crnt.s, B_start, B_end - 1, A_end);
        reinsert(solut_crnt.s, A_start, A_end - 1, B_end);
    } else {
        reinsert(solut_crnt.s, A_start, A_end - 1, B_end);
        reinsert(solut_crnt.s, B_start, B_end - 1, A_end);
    }
}

void solut_init(tSolution &solut, const tInfo &info) {
    solut.s.assign(static_cast<std::size_t>(info.dimen + 1), 0);
    solut.seq_dim = info.dimen + 1;
    solut.seq.assign(static_cast<std::size_t>(info.dimen + 1) * static_cast<std::size_t>(info.dimen + 1), tSeqInfo());
    solut.s_size = info.dimen + 1;
}

tSolution GILS_RVND(int Imax, int Iils, const std::vector<typeReal> &R, tInfo &info, bool verbose) {
    //! solutions
    tSolution sol_best;
    tSolution sol_partial;
    tSolution sol_crnt;

    //! aux variables
    int i;
    int index_;
    typeReal alpha;
    int iterILS;
#ifdef DEBUG
    bool fsb;
#endif

    int it;
    int Tit;

    solut_init(sol_best, info);
    solut_init(sol_partial, info);
    solut_init(sol_crnt, info);
    sol_best.cost = info.fmax;

    Tit = 0;
    for (i = 1; i <= Imax; i++) {

        index_ = info.rnd_at(info.rnd_index) + 1;
        info.rnd_index = info.rnd_index + 1;

        alpha = R[static_cast<std::size_t>(index_ - 1)];
        if (verbose) {
            std::cout << " [+] Local Search " << i << std::endl;
        }
        sol_crnt.s = construction(alpha, info);

#ifdef DEBUG
        // feasibility check
        fsb = false;
        is_feasible(sol_crnt, fsb);
        ASSERT(fsb, true);
#endif

        subseq_load(sol_crnt, info);
        sol_partial = sol_crnt;
        if (verbose) {
            std::cout << "         [+] Constructing Inital Solution.. " << sol_crnt.cost << std::endl;
            std::cout << "         [+] Looking for the best Neighbor.." << std::endl;
        }
        iterILS = 0;
        it = 0;

        while (iterILS < Iils) {
            RVND(sol_crnt, info, it);

            if (sol_crnt.cost < sol_partial.cost - std::numeric_limits<typeReal>::epsilon()) {
                sol_partial.s = sol_crnt.s;
                sol_partial.cost = sol_crnt.cost;
                iterILS = 0;
            }

            perturb(sol_crnt, sol_partial, info);
            subseq_load(sol_crnt, info);
#ifdef DEBUG
            fsb = false;
            is_feasible(sol_crnt, fsb);
            ASSERT(fsb, true);
#endif

            iterILS = iterILS + 1;
        }
        Tit = Tit + it;

        if (sol_partial.cost < sol_best.cost) {
            sol_best = sol_partial;
        }

        if (verbose) {
            std::cout << "         Current best solution cost: " << sol_best.cost << std::endl;
        }

    }

    (void)Tit;

    if (verbose) {
        for (std::size_t k = 0; k < sol_best.s.size(); k++) {
            std::cout << " " << sol_best.s[k];
        }
        std::cout << std::endl;
    }

    return sol_best;
}

int main(int argc, char **argv) {
    tInfo info;
    tSolution sol;
    std::vector<typeReal> R(26);
    int Iils;
    int Imax;
    int i;

    bool verbose;

    verbose = false;
    for (i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    data::load_matrix(info.cost, info.rnd);

    info.dimen = info.cost.n;

    {
        const float R_init[26] = {
            0.0f, 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f, 0.07f, 0.08f, 0.09f, 0.10f, 0.11f, 0.12f, 0.13f,
            0.14f, 0.15f, 0.16f, 0.17f, 0.18f, 0.19f, 0.20f, 0.21f, 0.22f, 0.23f, 0.24f, 0.25f
        };
        for (i = 0; i < 26; i++) {
            R[static_cast<std::size_t>(i)] = static_cast<typeReal>(R_init[i]);
        }
    }

    Iils = (100 < info.dimen ? 100 : info.dimen);
    Imax = 10;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::clock_t cpu_begin = std::clock();
    sol = GILS_RVND(Imax, Iils, R, info, verbose);
    std::clock_t cpu_end = std::clock();
    std::chrono::steady_clock::time_point end_ = std::chrono::steady_clock::now();

    std::printf(" COST:  %.10f\n", static_cast<double>(sol.cost));
    std::printf(" TIME:  %.7f\n", static_cast<double>(cpu_end - cpu_begin) / static_cast<double>(CLOCKS_PER_SEC));
    std::printf(" wall clock (s):  %.7f\n",
                std::chrono::duration_cast<std::chrono::duration<double> >(end_ - begin).count());

    return 0;
}
