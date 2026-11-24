#include "Data.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <string>
#include <vector>

// Indices of the subsequence info triple, as in main.lua's `info`.
static const int T = 1;
static const int C = 2;
static const int W = 3;

static const double EPSILON = 1e-15;

struct tInfo {
    std::vector<std::vector<double>> c;
    int dimension;
    std::vector<int> rnd;
    int rnd_index;
};

struct tSolution {
    std::vector<int> s;    // 1-based, s[0] unused
    std::vector<double> seq; // 1-based, seq[0] unused
    double cost;
};

// seq is a flat array indexed exactly like in the Lua version.
#define SEQ(i, j, k) seq[(3 * ((s_size) * ((i) - 1) + ((j) - 1)) + (k))]

static inline int sz(const std::vector<int> &v) { return (int)v.size() - 1; }

// ---------------------------------------------------------------------------
// Lua-compatible number formatting (%.14g), so the printed output matches
// main.lua's.
// ---------------------------------------------------------------------------
static std::string num_str(double v) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.14g", v);
    return std::string(buf);
}

// ---------------------------------------------------------------------------
// math.random(m, n) equivalent
// ---------------------------------------------------------------------------
static inline int rand_range(int m, int n) {
    if (n < m) return m;
    return m + (int)(rand() % (unsigned)(n - m + 1));
}

void s_print(const tSolution &solut) {
    for (int i = 1; i <= sz(solut.s); i++) {
        printf("%d ", solut.s[i]);
    }
    printf("\n");
}

void table_print(const std::vector<int> &tbl) {
    for (int i = 1; i <= sz(tbl); i++) {
        printf("%d ", tbl[i]);
    }
    printf("\n");
}

void matrix_print(const tInfo &info) {
    for (int i = 1; i < (int)info.c.size(); i++) {
        for (int j = 1; j < (int)info.c.size(); j++) {
            printf("%s ", num_str(info.c[i][j]).c_str());
        }
        printf("\n");
    }
}

void seq_print(const tSolution &solut, const tInfo &info) {
    const std::vector<double> &seq = solut.seq;
    const int s_size = info.dimension + 1;

    for (int i = 1; i <= info.dimension + 1; i++) {
        for (int j = i; j <= info.dimension + 1; j++) {
            for (int k = 1; k <= 3; k++) {
                printf("%s ", num_str(SEQ(i, j, k)).c_str());
            }
            printf("| ");
        }
        printf("\n");
    }
}

void solut_clone(const tSolution &src, tSolution &dest) {
    dest.s.resize(src.s.size());
    for (int i = 1; i <= sz(src.s); i++) {
        dest.s[i] = src.s[i];
    }
    dest.cost = src.cost;
}

void subseq_fill(std::vector<double> &seq, const tInfo &info) {
    const int s_size = info.dimension + 1;
    for (int i = 1; i <= info.dimension + 1; i++) {
        for (int j = i; j <= info.dimension + 1; j++) {
            for (int k = 1; k <= 3; k++) {
                SEQ(i, j, k) = 0.0;
            }
        }
    }
}

void subseq_load(tSolution &solut, const tInfo &info) {
    const int dimen = info.dimension;
    std::vector<double> &seq = solut.seq;

    const std::vector<std::vector<double>> &c = info.c;
    const std::vector<int> &s = solut.s;
    const int s_size = dimen + 1;

    for (int i = 1; i <= s_size; i++) {
        const int k = 1 - i - (i != 1 ? 0 : 1);

        SEQ(i, i, T) = 0.0;
        SEQ(i, i, C) = 0.0;
        SEQ(i, i, W) = (i != 1 ? 1 : 0);

        double some_T = 0.0;
        double some_C = 0.0;

        int s_j_prev = s[i];

        for (int j = i + 1; j <= s_size; j++) {
            const int s_j = s[j];

            // set T
            const double cost_T = c[s_j_prev][s_j] + some_T;
            SEQ(i, j, T) = cost_T;

            // set C
            const double cost_C = cost_T + some_C;
            SEQ(i, j, C) = cost_C;

            // set W
            SEQ(i, j, W) = j + k;

            s_j_prev = s_j;
            some_T = cost_T;
            some_C = cost_C;
        }
    }

    solut.cost = SEQ(1, dimen + 1, C) - EPSILON;
}

int partition(std::vector<int> &arr, int left, int right, int r, const tInfo &info) {
    const int pivotValue = arr[right];
    int i = left - 1;

    for (int j = left; j <= right - 1; j++) {
        if (info.c[r][arr[j]] < info.c[r][pivotValue]) {
            i = i + 1;
            std::swap(arr[i], arr[j]);
        }
    }

    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

void quicksort(std::vector<int> &arr, int left, int right, int r, const tInfo &info) {
    if (left < right) {
        const int pivotIndex = partition(arr, left, right, r, info);
        quicksort(arr, left, pivotIndex - 1, r, info);
        quicksort(arr, pivotIndex + 1, right, r, info);
    }
}

void sort(std::vector<int> &arr, int r, const tInfo &info) {
    quicksort(arr, 1, sz(arr), r, info);
}

// table.insert(t, pos, e) / table.remove(t, pos) on a 1-based vector
static inline void tbl_insert(std::vector<int> &t, int pos, int e) {
    t.insert(t.begin() + pos, e);
}

static inline int tbl_remove(std::vector<int> &t, int pos) {
    const int e = t[pos];
    t.erase(t.begin() + pos);
    return e;
}

std::vector<int> construction(double alpha, tInfo &info) {
    std::vector<int> s(1, 0);
    s.push_back(1);

    std::vector<int> cList(1, 0);
    for (int i = 2; i <= info.dimension; i++) {
        cList.push_back(i);
    }

    int r = 1;

    while (sz(cList) > 0) {
        sort(cList, r, info);

        const int range = (int)floor(sz(cList) * alpha) + 1;
        int i = rand_range(1, range);
        i = info.rnd[info.rnd_index] + 1;
        info.rnd_index = info.rnd_index + 1;

        const int c = tbl_remove(cList, i);
        s.push_back(c);
        r = c;
    }
    s.push_back(1);

    return s;
}

void swap(std::vector<int> &s, int i, int j) {
    std::swap(s[i], s[j]);
}

void reverse(std::vector<int> &s, int i, int j) {
    int l = j;
    for (int k = i; k <= (int)floor((j + i) / 2.0); k++) {
        swap(s, k, l);
        l = l - 1;
    }
}

void reinsert(std::vector<int> &s, int i, int j, int pos) {
    if (i < pos) {
        for (int k = i; k <= j; k++) {
            const int e = s[i];
            tbl_insert(s, pos, e);
            tbl_remove(s, i);
        }
    } else {
        for (int k = i; k <= j; k++) {
            const int e = s[j];
            tbl_remove(s, j);
            tbl_insert(s, pos, e);
        }
    }
}

bool search_swap(tSolution &solut, const tInfo &info) {
    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;
    const int dimen = info.dimension;

    double cost_concat_1 = 0.0;
    double cost_concat_2 = 0.0;
    double cost_concat_3 = 0.0;
    double cost_concat_4 = 0.0;
    double cost_new = 0.0;

    const std::vector<int> &s = solut.s;
    const std::vector<double> &seq = solut.seq;
    const std::vector<std::vector<double>> &c = info.c;

    const int s_size = dimen + 1;

    const double EP = EPSILON;

    for (int i = 2; i <= dimen - 1; i++) {
        const int i_prev = i - 1;
        const int i_next = i + 1;

        const int s_i = s[i];
        const int s_i_next = s[i_next];

        cost_concat_1 =                 SEQ(1, i_prev, T) + c[s[i_prev]][s_i_next];
        cost_concat_2 = cost_concat_1 + SEQ(i, i_next, T) + c[s_i][     s[i_next + 1]];

        cost_new = SEQ(1, i_prev, C)                                                    +           //       1st subseq
        SEQ(i, i_next, W)               * (cost_concat_1) + c[s_i_next][s_i]  +           // concat 2nd subseq
        SEQ(i_next + 1, s_size, W)   * (cost_concat_2) + SEQ(i_next + 1, s_size, C);  // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - EP;
            I = i;
            J = i_next;
        }

        for (int j = i_next + 1; j <= dimen; j++) {
            const int j_prev = j - 1;
            const int j_next = j + 1;
            const int s_j = s[j];

            cost_concat_1 =                 SEQ(1, i_prev, T)       + c[s[i_prev]][s_j];
            cost_concat_2 = cost_concat_1                     + c[s_j][s_i_next];
            cost_concat_3 = cost_concat_2 + SEQ(i_next, j_prev, T)  + c[s[j_prev]][s_i];
            cost_concat_4 = cost_concat_3                           + c[s_i][s[j_next]];

            cost_new = SEQ(1, i_prev, C)                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            SEQ(i_next, j_prev, W)      * cost_concat_2 + SEQ(i_next, j_prev, C) +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            SEQ(j_next, s_size, W) * cost_concat_4 + SEQ(j_next, s_size, C);   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - EP;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < SEQ(1, s_size, C) - EPSILON) {
        swap(solut.s, I, J);
        subseq_load(solut, info);
        return true;
    }

    return false;
}

bool search_two_opt(tSolution &solut, const tInfo &info) {
    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;
    const int dimen = info.dimension;

    double cost_concat_1 = 0.0;
    double cost_concat_2 = 0.0;
    double cost_new = 0.0;

    const std::vector<int> &s = solut.s;
    const std::vector<double> &seq = solut.seq;
    const std::vector<std::vector<double>> &c = info.c;

    const double EP = EPSILON;

    const int s_size = dimen + 1;

    for (int i = 2; i <= dimen - 1; i++) {
        const int i_prev = i - 1;
        double rev_seq_cost = SEQ(i, i + 1, T);

        int s_j_prev = s[i + 1];

        for (int j = i + 2; j <= dimen; j++) {
            const int j_next = j + 1;
            const int s_j = s[j];

            rev_seq_cost = rev_seq_cost + c[s_j_prev][s_j] * (SEQ(i, j, W) - 1.0);

            cost_concat_1 =                 SEQ(1, i_prev, T)   + c[s_j][     s[i_prev]];
            cost_concat_2 = cost_concat_1 + SEQ(i, j, T)        + c[s[j_next]][s[i]];

            cost_new = SEQ(1, i_prev, C)                                                        +   //   1st subseq
                    SEQ(i, j, W)                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
                    SEQ(j_next, s_size, W) * cost_concat_2 + SEQ(j_next, s_size, C);    // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - EP;
                I = i;
                J = j;
            }

            s_j_prev = s_j;
        }
    }

    if (cost_best < SEQ(1, s_size, C) - EP) {
        reverse(solut.s, I, J);
        subseq_load(solut, info);
        return true;
    }

    return false;
}

bool search_reinsertion(tSolution &solut, const tInfo &info, int opt) {
    double cost_best = std::numeric_limits<double>::infinity();
    int I = -1;
    int J = -1;
    int POS = -1;
    const int dimen = info.dimension;

    double cost_concat_1 = 0.0;
    double cost_concat_2 = 0.0;
    double cost_concat_3 = 0.0;
    double cost_new = 0.0;

    const std::vector<int> &s = solut.s;
    const std::vector<double> &seq = solut.seq;
    const std::vector<std::vector<double>> &c = info.c;
    const int s_size = dimen + 1;

    const double EP = EPSILON;

    for (int i = 2; i <= dimen - opt + 1; i++) {
        const int j = opt + i - 1;
        const int i_prev = i - 1;
        const int j_next = j + 1;

        const int s_i = s[i];
        const int s_j = s[j];

        const int s_i_prev = s[i_prev];
        const int s_j_next = s[j_next];

        for (int k = 1; k <= i_prev - 1; k++) {
            const int k_next = k + 1;

            cost_concat_1 =                 SEQ(1, k, T)            + c[s[k]][s_i];
            cost_concat_2 = cost_concat_1 + SEQ(i, j, T)            + c[s_j][s[k_next]];
            cost_concat_3 = cost_concat_2 + SEQ(k_next, i_prev, T)  + c[s_i_prev][s_j_next];

            cost_new = SEQ(1, k, C)                                                             +   //       1st subseq
            SEQ(i, j, W)                * cost_concat_1 + SEQ(i, j, C)                  +   // concat 2nd subseq (reinserted seq)
            SEQ(k_next, i_prev, W)      * cost_concat_2 + SEQ(k_next, i_prev, C)        +   // concat 3rd subseq
            SEQ(j_next, s_size, W) * cost_concat_3     + SEQ(j_next, s_size, C);       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - EP;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (int k = i + opt; k <= dimen; k++) {
            const int k_next = k + 1;

            cost_concat_1 =                 SEQ(1, i_prev, T)   + c[s_i_prev][s_j_next];
            cost_concat_2 = cost_concat_1 + SEQ(j_next, k, T)   + c[s[k]][     s_i];
            cost_concat_3 = cost_concat_2 + SEQ(i, j, T)        + c[s_j][     s[k_next]];

            cost_new = SEQ(1, i_prev, C)                                                        +   //       1st subseq
                    SEQ(j_next, k, W)           * cost_concat_1 + SEQ(j_next, k, C)             +   // concat 2nd subseq
                    SEQ(i, j, W)                * cost_concat_2 + SEQ(i, j, C)                  +   // concat 3rd subseq (reinserted seq)
                    SEQ(k_next, s_size, W) * cost_concat_3     + SEQ(k_next, s_size, C);       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - EP;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost) {
        reinsert(solut.s, I, J, POS + 1);
        subseq_load(solut, info);

        if (cost_best != solut.cost) {
            printf("ERROR\n");
            exit(1);
        }
        return true;
    }

    return false;
}

void RVND(tSolution &solut, tInfo &info) {
    const int SWAP        = 0;
    const int REINSERTION = 1;
    const int OR_OPT_2    = 2;
    const int OR_OPT_3    = 3;
    const int TWO_OPT     = 4;

    const int init_list[] = {0, SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    std::vector<int> neighbd_list(init_list, init_list + 6);

    while (sz(neighbd_list) > 0) {
        int index = rand_range(1, sz(neighbd_list));

        index = info.rnd[info.rnd_index] + 1;
        info.rnd_index = info.rnd_index + 1;

        const int neighbd = neighbd_list[index];

        bool improve = false;

        if (neighbd == SWAP) {
            improve = search_swap(solut, info);
        } else if (neighbd == REINSERTION) {
            improve = search_reinsertion(solut, info, REINSERTION);
        } else if (neighbd == OR_OPT_2) {
            improve = search_reinsertion(solut, info, OR_OPT_2);
        } else if (neighbd == OR_OPT_3) {
            improve = search_reinsertion(solut, info, OR_OPT_3);
        } else if (neighbd == TWO_OPT) {
            improve = search_two_opt(solut, info);
        }

        if (improve == true) {
            neighbd_list.assign(init_list, init_list + 6);
        } else {
            tbl_remove(neighbd_list, index);
        }
    }
}

std::vector<int> perturb(const tSolution &solut, tInfo &info) {
    std::vector<int> s(solut.s);

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = (int)floor(sz(s) / 10.0);
    size_max = size_max >= 2 ? size_max : 2;
    const int size_min = 2;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        A_start = rand_range(2, sz(s) - 1 - size_max);
        A_end = A_start + rand_range(size_min, size_max);

        B_start = rand_range(2, sz(s) - 1 - size_max);
        B_end = B_start + rand_range(size_min, size_max);

        A_start = info.rnd[info.rnd_index] + 1;
        info.rnd_index = info.rnd_index + 1;
        A_end = A_start + info.rnd[info.rnd_index];
        info.rnd_index = info.rnd_index + 1;

        B_start = info.rnd[info.rnd_index] + 1;
        info.rnd_index = info.rnd_index + 1;
        B_end = B_start + info.rnd[info.rnd_index];
        info.rnd_index = info.rnd_index + 1;
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

void GILS_RVND(int Imax, int Iils, const std::vector<double> &R, tInfo &info, bool verbose) {

    const int seq_n = 3 * (info.dimension + 1) * (info.dimension + 1);

    tSolution solut_partial;
    solut_partial.s.assign(1, 0);
    solut_partial.seq.assign(seq_n + 1, 0.0);

    tSolution solut_crnt;
    solut_crnt.s.assign(1, 0);
    solut_crnt.seq.assign(seq_n + 1, 0.0);

    tSolution solut_best;
    solut_best.s.assign(1, 0);
    solut_best.seq.assign(seq_n + 1, 0.0);

    subseq_fill(solut_partial.seq, info);
    solut_partial.cost = 0;

    subseq_fill(solut_crnt.seq, info);
    solut_crnt.cost = 0;

    subseq_fill(solut_best.seq, info);
    solut_best.cost = 0;

    solut_best.cost = std::numeric_limits<double>::infinity();

    for (int i = 1; i <= Imax; i++) {
        const int Rsz = (int)R.size() - 1;
        double alpha = R[rand_range(1, Rsz)];
        alpha = R[info.rnd[info.rnd_index] + 1];
        info.rnd_index = info.rnd_index + 1;

        if (verbose) {
            printf("[+] Local Search\t%d\n", i);
            printf("\t[+] Constructing Inital Solution..\n");
        }
        solut_crnt.s = construction(alpha, info);
        subseq_load(solut_crnt, info);
        if (verbose) {
            s_print(solut_crnt);
            printf("\tConstruction cost  \t%s\n", num_str(solut_crnt.cost).c_str());
        }

        solut_clone(solut_crnt, solut_partial);

        if (verbose) {
            printf("\t[+] Looking for the best Neighbor..\n");
        }
        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(solut_crnt, info);
            if (solut_crnt.cost < solut_partial.cost - EPSILON) {
                solut_partial.cost = solut_crnt.cost - EPSILON;
                solut_partial.s = solut_crnt.s;
                iterILS = 0;
            }

            solut_crnt.s = perturb(solut_partial, info);
            subseq_load(solut_crnt, info);

            iterILS = iterILS + 1;
        }

        subseq_load(solut_partial, info);

        if (solut_partial.cost < solut_best.cost) {
            solut_clone(solut_partial, solut_best);
        }
        if (verbose) {
            printf("\tCurrent best solution cost: \t%s\n", num_str(solut_best.cost).c_str());
        }
    }

    printf("COST: \t%s\n", num_str(solut_best.cost).c_str());
}

int main(int argc, char **argv) {
    tInfo info;
    info.rnd.assign(1, 0);
    info.rnd_index = 1;

    int a = 0;
    info.dimension = readData(info.c, info.rnd, a);
    printf("%d\n", info.rnd[a]);
    srand((unsigned)time(NULL));

    const int Imax = 10;
    const int Iils = info.dimension < 100 ? info.dimension : 100;
    const double R_vals[] = {0.0,
               0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
               0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    const std::vector<double> R(R_vals, R_vals + sizeof(R_vals) / sizeof(R_vals[0]));

    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    const clock_t start = clock();
    GILS_RVND(Imax, Iils, R, info, verbose);

    printf("TIME: \t%s\n", num_str((double)(clock() - start) / CLOCKS_PER_SEC).c_str());

    return 0;
}
