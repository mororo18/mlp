#include "Data.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>

struct Info {
    std::vector<std::vector<double> > c;
    std::vector<int> rnd;
    int dimension;
    int T;
    int C;
    int W;
    int rnd_index;
};

static const double EPSILON = std::numeric_limits<double>::epsilon();
static const double MAX_VALUE = std::numeric_limits<double>::max();

/* prints a double the way console.log does */
static std::string num_to_string(double value) {
    char buffer[64];

    if (value == std::floor(value) && std::fabs(value) < 1e21) {
        std::snprintf(buffer, sizeof(buffer), "%.0f", value);
        return std::string(buffer);
    }

    for (int precision = 1; precision <= 17; precision++) {
        std::snprintf(buffer, sizeof(buffer), "%.*g", precision, value);
        if (std::atof(buffer) == value) {
            break;
        }
    }
    return std::string(buffer);
}

static int to_1D(int x, int y, int z, int n) {
    return (x*n + y)*3 + z;
}

static std::vector<double> subseq_fill(int dimension) {
    std::vector<double> seq((dimension+1)*(dimension+1)*3);
    return seq;
}

static void subseq_load(const std::vector<int> & s, std::vector<double> & seq, Info & info) {

    for (int i = 0; i < info.dimension+1; i++) {
        int k = 1 - i - (i != 0 ? 0 : 1);

        seq[to_1D(i, i, info.T, info.dimension)] = 0.0;
        seq[to_1D(i, i, info.C, info.dimension)] = 0.0;
        seq[to_1D(i, i, info.W, info.dimension)] = (i != 0 ? 1.0 : 0.0);

        for (int j = i+1; j < info.dimension+1; j++) {
            int j_prev = j-1;
            seq[to_1D(i, j, info.T, info.dimension)] = info.c[s[j_prev]][s[j]] + seq[to_1D(i, j_prev, info.T, info.dimension)];
            seq[to_1D(i, j, info.C, info.dimension)] = seq[to_1D(i, j, info.T, info.dimension)] + seq[to_1D(i, j_prev, info.C, info.dimension)];
            seq[to_1D(i, j, info.W, info.dimension)] = j + k;
        }
    }
}

static int partition(std::vector<int> & arr, int left, int right, int r, Info & info) {
    int pivotValue = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (info.c[r][arr[j]] < info.c[r][pivotValue]) {
            i++;
            std::swap(arr[i], arr[j]); // Swap elements
        }
    }
    std::swap(arr[i + 1], arr[right]); // Swap pivot
    return i + 1;
}

static void quicksort(std::vector<int> & arr, int left, int right, int r, Info & info) {
    if (left < right) {
        int pivotIndex = partition(arr, left, right, r, info);
        quicksort(arr, left, pivotIndex - 1, r, info);
        quicksort(arr, pivotIndex + 1, right, r, info);
    }
}

static void sort(std::vector<int> & arr, int r, Info & info) {
    quicksort(arr, 0, (int)arr.size() - 1, r, info);
}

static std::vector<int> construction(double alpha, Info & info) {
    (void)alpha;
    std::vector<int> s(1, 0);
    std::vector<int> cList(info.dimension-1);
    for (int i = 0; i < info.dimension-1; i++) {
        cList[i] = i + 1;
    }

    int r = 0;
    while (cList.size() > 0) {
        sort(cList, r, info);

        int i = info.rnd[info.rnd_index++];
        int c = cList[i];
        cList.erase(cList.begin() + i);
        s.push_back(c);
        r = c;
    }

    s.resize(info.dimension + 1);
    s[info.dimension] = 0;

    return s;
}

static void swap(std::vector<int> & s, int i, int j) {
    std::swap(s[i], s[j]);
}

static void reverse(std::vector<int> & s, int i, int j) {
    for (int f = i, b = j; f <= (i+j)/2; f++, b--) {
        swap(s, f, b);
    }
}

static void reinsert(std::vector<int> & s, int i, int j, int pos) {
    int sz = j - i + 1;
    std::vector<int> sub(s.begin() + i, s.begin() + i + sz);
    if (i < pos) {
        s.erase(s.begin() + i, s.begin() + i + sz);
        s.insert(s.begin() + (pos - sz), sub.begin(), sub.end());
    } else {
        s.erase(s.begin() + i, s.begin() + i + sz);
        s.insert(s.begin() + pos, sub.begin(), sub.end());
    }
}

static bool search_swap(std::vector<int> & s, std::vector<double> & seq, Info & info) {
    double cost_best = MAX_VALUE;
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    double cost_concat_4;
    int I = -1;
    int J = -1;

    for (int i = 1; i < info.dimension-1; i++) {
        int i_prev = i - 1;
        int i_next = i + 1;

        // immediate nodes case

        cost_concat_1 =                 seq[to_1D(0, i_prev, info.T, info.dimension)] + info.c[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[to_1D(i, i_next, info.T, info.dimension)] + info.c[s[i]][s[i_next+1]];

        cost_new = seq[to_1D(0, i_prev, info.C, info.dimension)]                                                            // 1st subseq
            + seq[to_1D(i, i_next, info.W, info.dimension)]           * (cost_concat_1) + info.c[s[i_next]][s[i]]              // concat 2nd subseq
            + seq[to_1D(i_next+1, info.dimension, info.W, info.dimension)] * (cost_concat_2) + seq[to_1D(i_next+1, info.dimension, info.C, info.dimension)];  // concat 3rd subseq

        if(cost_new < cost_best){
            cost_best = cost_new - EPSILON;
            I = i;
            J = i_next;
        }

        for(int j = i_next+1; j < info.dimension; j++){
            int j_next = j+1;
            int j_prev = j-1;

            cost_concat_1 = seq[to_1D(0, i_prev, info.T, info.dimension)] + info.c[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1 + info.c[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[to_1D(i_next, j_prev, info.T, info.dimension)] + info.c[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3 + info.c[s[i]][s[j_next]];

            cost_new = seq[to_1D(0, i_prev, info.C, info.dimension)]                                                        /* first subseq */
                + cost_concat_1                                                             /* concatenate second subseq (single node) */
                + seq[to_1D(i_next, j_prev, info.W, info.dimension)] * cost_concat_2 + seq[to_1D(i_next, j_prev, info.C, info.dimension)]           /* concatenate third subseq */
                + cost_concat_3                                                             /* concatenate fourth subseq (single node) */
                + seq[to_1D(j_next, info.dimension, info.W, info.dimension)] * cost_concat_4 + seq[to_1D(j_next, info.dimension, info.C, info.dimension)];    /* concatenate fifth subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
            }

        }

    }


    if(cost_best < seq[to_1D(0, info.dimension, info.C, info.dimension)] - EPSILON){
        swap(s, I, J);
        subseq_load(s, seq, info);
        return true;
    }

    return false;
}

static bool search_two_opt(std::vector<int> & s, std::vector<double> & seq, Info & info) {
    int I = -1;
    int J = -1;
    double cost_best = MAX_VALUE;
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;

    for(int i = 1; i < info.dimension-1; i++){
        int i_prev = i -1;
        double rev_seq_cost = seq[to_1D(i, i+1, info.T, info.dimension)];

        for(int j = i+2; j < info.dimension; j++){
            int j_next = j+1;
            int j_prev = j-1;

            rev_seq_cost += info.c[s[j_prev]][s[j]] * (seq[to_1D(i, j, info.W, info.dimension)]-1);

            cost_concat_1 =                 seq[to_1D(0, i_prev, info.T, info.dimension)] + info.c[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[to_1D(i, j, info.T, info.dimension)] + info.c[s[j_next]][s[i]];

            cost_new = seq[to_1D(0, i_prev, info.C, info.dimension)]                                                        /*        1st subseq */
                + seq[to_1D(i, j, info.W, info.dimension)]              * cost_concat_1 + rev_seq_cost                  /* concat 2nd subseq (reversed seq) */
                + seq[to_1D(j_next, info.dimension, info.W, info.dimension)] * cost_concat_2 + seq[to_1D(j_next, info.dimension, info.C, info.dimension)];    /* concat 3rd subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if(cost_best < seq[to_1D(0, info.dimension, info.C, info.dimension)] - EPSILON){
        reverse(s, I, J);
        subseq_load(s, seq, info);
        return true;
    }
    return false;
}

static bool search_reinsertion(std::vector<int> & s, std::vector<double> & seq, Info & info, int opt) {
    double cost_best = MAX_VALUE;
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    int I = -1;
    int J = -1;
    int POS = -1;

    for (int i = 1; i < info.dimension - opt + 1; i++) {
        int j = opt + i - 1;
        int i_prev = i-1;
        int j_next = j+1;

        // k -> reinsertion places
        for (int k = 0; k < i_prev; k++) {
            int k_next = k+1;

            cost_concat_1 = seq[to_1D(0, k, info.T, info.dimension)] + info.c[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[to_1D(i, j, info.T, info.dimension)] + info.c[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[to_1D(k_next, i_prev, info.T, info.dimension)] + info.c[s[i_prev]][s[j_next]];

            cost_new = seq[to_1D(0, k, info.C, info.dimension)]                                                             /*        1st subseq */
                + seq[to_1D(i, j, info.W, info.dimension)]              * cost_concat_1 + seq[to_1D(i, j, info.C, info.dimension)]                  /* concat 2nd subseq (reinserted seq) */
                + seq[to_1D(k_next, i_prev, info.W, info.dimension)]    * cost_concat_2 + seq[to_1D(k_next, i_prev, info.C, info.dimension)]        /* concat 3rd subseq */
                + seq[to_1D(j_next, info.dimension, info.W, info.dimension)] * cost_concat_3 + seq[to_1D(j_next, info.dimension, info.C, info.dimension)];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (int k = i+opt; k < info.dimension; k++) {
            int k_next = k+1;

            cost_concat_1 = seq[to_1D(0, i_prev, info.T, info.dimension)] + info.c[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[to_1D(j_next, k, info.T, info.dimension)] + info.c[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[to_1D(i, j, info.T, info.dimension)] + info.c[s[j]][s[k_next]];

            cost_new = seq[to_1D(0, i_prev, info.C, info.dimension)]                                                        /*      1st subseq */
                + seq[to_1D(j_next, k, info.W, info.dimension)]         * cost_concat_1 + seq[to_1D(j_next, k, info.C, info.dimension)]             /* concat 2nd subseq */
                + seq[to_1D(i, j, info.W, info.dimension)]              * cost_concat_2 + seq[to_1D(i, j, info.C, info.dimension)]                  /* concat 3rd subseq (reinserted seq) */
                + seq[to_1D(k_next, info.dimension, info.W, info.dimension)] * cost_concat_3 + seq[to_1D(k_next, info.dimension, info.C, info.dimension)];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if(cost_best < seq[to_1D(0, info.dimension, info.C, info.dimension)] - EPSILON){
        reinsert(s, I, J, POS+1);
        subseq_load(s, seq, info);

        return true;
    }

    return false;
}

static void RVND(std::vector<int> & s, std::vector<double> & subseq, Info & info) {
    const int SWAP        = 0;
    const int REINSERTION = 1;
    const int OR_OPT_2    = 2;
    const int OR_OPT_3    = 3;
    const int TWO_OPT     = 4;

    const int initial[5] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    std::vector<int> neighbd_list(initial, initial + 5);

    bool improve = false;

    while (neighbd_list.size() > 0) {
        int i = info.rnd[info.rnd_index++];
        int neighbd = neighbd_list[i];

        switch (neighbd) {
            case SWAP:
                improve = search_swap(s, subseq, info);
                break;
            case REINSERTION:
                improve = search_reinsertion(s, subseq, info, REINSERTION);
                break;
            case OR_OPT_2:
                improve = search_reinsertion(s, subseq, info, OR_OPT_2);
                break;
            case OR_OPT_3:
                improve = search_reinsertion(s, subseq, info, OR_OPT_3);
                break;
            case TWO_OPT:
                improve = search_two_opt(s, subseq, info);
                break;
        }

        if (improve) {
            neighbd_list.assign(initial, initial + 5);
        } else {
            neighbd_list.erase(neighbd_list.begin() + i);
        }
    }

}

static std::vector<int> perturb(const std::vector<int> & sl, Info & info) {
    std::vector<int> s(sl);

    int A_start = 1, A_end = 1;
    int B_start = 1, B_end = 1;

    int size_max = (int)std::floor(sl.size()/10.0);
    size_max = (size_max >= 2 ? size_max : 2);
    int size_min = 2;
    (void)size_max;
    (void)size_min;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = info.rnd[info.rnd_index++];
        A_end = A_start + info.rnd[info.rnd_index++];

        B_start = info.rnd[info.rnd_index++];
        B_end = B_start + info.rnd[info.rnd_index++];
    }

    if(A_start < B_start){
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    }else{
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    return s;

}

static void print_solution(const std::vector<int> & s) {
    std::cout << "[";
    for (size_t i = 0; i < s.size(); i++) {
        std::cout << (i == 0 ? " " : ", ") << s[i];
    }
    std::cout << (s.size() > 0 ? " ]" : "]") << std::endl;
}

static void GILS_RVND(int Iils, int Imax, const std::vector<double> & R, Info & info, bool verbose) {

    std::vector<double> subseq = subseq_fill(info.dimension);
    std::vector<int> s_best;
    double cost_best = MAX_VALUE;

    for (int i = 0; i < Imax; i++) {
        double alpha = R[info.rnd[info.rnd_index++]];
        if (verbose) {
            std::cout << "[+] Local Search  " << i+1 << std::endl;
            std::cout << "\t[+] Constructing Inital Solution.." << std::endl;
        }
        std::vector<int> s = construction(alpha, info);
        std::vector<int> sl(s);

        subseq_load(s, subseq, info);
        double rvnd_cost_best = subseq[to_1D(0, info.dimension, info.C, info.dimension)] - EPSILON;
        if (verbose) {
            std::cout << "Construction cost " << num_to_string(rvnd_cost_best) << std::endl;
        }
        int iterILS = 0;

        if (verbose) {
            std::cout << "\t[+] Looking for the best Neighbor.." << std::endl;
        }
        while (iterILS < Iils) {
            RVND(s, subseq, info);
            double rvnd_cost_crnt = subseq[to_1D(0, info.dimension, info.C, info.dimension)] - EPSILON;
            if (rvnd_cost_crnt < rvnd_cost_best) {
                rvnd_cost_best = rvnd_cost_crnt;
                sl = s;
                iterILS = 0;
            }

            s = perturb(sl, info);
            subseq_load(s, subseq, info);
            iterILS++;
        }


        subseq_load(sl, subseq, info);
        double sl_cost = subseq[to_1D(0, info.dimension, info.C, info.dimension)] - EPSILON;

        if (sl_cost < cost_best) {
            s_best = sl;
            cost_best = sl_cost;
        }
        if (verbose) {
            std::cout << "\tcurrent best solution cost " << num_to_string(cost_best) << std::endl;
            std::cout << std::endl;
        }
    }

    print_solution(s_best);
    std::cout << "COST:  " << num_to_string(cost_best) << std::endl;
}

static double cpu_usage_seconds() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    double user = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec / 1e6;
    double system = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec / 1e6;
    return user + system;
}

static double wall_clock_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

int main(int argc, char ** argv) {
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    int dimension;
    std::vector<std::vector<double> > c;
    std::vector<int> rnd;

    InfoLoadResult t = info_load(c);
    dimension = t.dimension;
    rnd = t.rnd;
    int Iils = std::min(dimension, 100);
    const int Imax = 10;
    const double R_values[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    const std::vector<double> R(R_values, R_values + sizeof(R_values)/sizeof(R_values[0]));

    Info info;
    info.c = c;
    info.rnd = rnd;
    info.dimension = dimension;
    info.T = 0;
    info.C = 1;
    info.W = 2;
    info.rnd_index = 0;

    double cpuStart = cpu_usage_seconds();
    double start = wall_clock_seconds();
    GILS_RVND(Iils, Imax, R, info, verbose);
    double end = wall_clock_seconds();
    double cpuUsage = cpu_usage_seconds() - cpuStart;

    std::cout << "TIME:  " << num_to_string(cpuUsage) << std::endl;
    std::cout << "wall clock (s):  " << num_to_string(end - start) << std::endl;

    return 0;
}
