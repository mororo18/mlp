#include <cstdio>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <vector>
#include "data.hpp"
#include "types.hpp"

#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

Real R_table(int i){
    static const Real table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

void print_s(const std::vector<int> & s, int sz) {

    for (int i = 0; i < sz; i++)
        printf("%d ", s[i]);
    printf("\n");
}

int partition(std::vector<int> & arr, int left, int right, const tData & data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (data.cost[r][arr[j]] < data.cost[r][pivot]) {
            i++;
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    int temp = arr[i + 1];
    arr[i + 1] = arr[right];
    arr[right] = temp;
    return i + 1;
}

void quicksort(std::vector<int> & arr, int left, int right, const tData & data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

void sort(std::vector<int> & arr, int len, int r, const tData & data) {
    quicksort(arr, 0, len - 1, data, r);
}

#ifndef NDEBUG
bool feasible(const std::vector<int> & s, int sz) {
    assert(sz > 0);
    std::vector<int> count(sz, 0);

    for (int i = 0; i < sz; i++) {
        assert(s[i] >= 0 && s[i] < sz);
    }

    for (int i = 0; i < sz; i++)
        count[s[i]] += 1;

    for (int i = 0; i < sz; i++)
        if (count[i] != 1) {
            return false;
        }

    return true;
}
#endif

void construct(std::vector<int> & ret, const Real alpha, tData & data){

    std::vector<int> s(data.dimen+1, 0);
    std::vector<int> cL(data.dimen-1, 0);
    int cL_size = data.dimen-1;

    for (int i = 1; i < data.dimen; ++i) {
        cL[i-1] = i;
    }

    int r = 0;
    for (int j = 1; j < data.dimen; ++j) {
        sort(cL, cL_size, r, data);

        int index = data.rnd[data.rnd_index++];
        int c = cL[index];
        s[j] = c;
        r = c;
        std::copy(cL.begin()+index+1, cL.begin()+cL_size, cL.begin()+index);
        cL_size--;
    }

    std::copy(s.begin(), s.begin()+(data.dimen+1), ret.begin());
}

void swap(std::vector<int> & vec, int i, int j){
    int tmp = vec[i]; vec[i] = vec[j]; vec[j] = tmp;
}

void reverse(std::vector<int> & vec, int i, int j){
    int m = (i+j)/2;
    for (int first = i, last = j;
            first <= m;
            first++, last--) {

        swap(vec, first, last);
    }
}

void reinsert(std::vector<int> & vec, int i, int j, int pos){
    std::vector<int> seq(vec.begin()+i, vec.begin()+j+1);

    if(pos < i){
        int sz = i-pos;
        std::memmove(vec.data() + j+1-sz, vec.data() + pos, sizeof(int)*sz);
        std::copy(seq.begin(), seq.end(), vec.begin()+pos);
    }else{
        int sz = pos-j-1;
        std::memmove(vec.data() + i, vec.data() + j+1, sizeof(int)*sz);
        std::copy(seq.begin(), seq.end(), vec.begin()+i+sz);
    }
}

void update_subseq_info_matrix(tSolution & solut, const tData & data){
    int i, j, j_prev, k;
    int from = 0;
    bool t;
    for (i = 0; i < data.dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut.seq_set_T(i, i, 0.0);
        solut.seq_set_C(i, i, 0.0);
        solut.seq_set_W(i, i, (Real)!(i == 0));

        for (j = i+1; j < data.dimen+1; j++) {
            j_prev = j-1;

            Real T = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq_get_T(i, j_prev);
            solut.seq_set_T(i, j, T);

            Real C = solut.seq_get_T(i, j) + solut.seq_get_C(i, j_prev);
            solut.seq_set_C(i, j, C);

            Real W = j + k;
            solut.seq_set_W(i, j, W);

        }
        from += t;
    }

    solut.cost = solut.seq_get_C(0, data.dimen);
}

bool search_swap(tSolution & solut, const tData & data) {
    Real cost_new,
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    Real cost_best = DBL_MAX;
    int i, j, j_prev, j_next, i_prev, i_next;
    int I;
    int J;

    for (i = 1; i < data.dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut.seq_get_T(0, i_prev) + data.cost[solut.s[i_prev]][solut.s[i_next]];
        cost_concat_2 = cost_concat_1 + solut.seq_get_T(i, i_next)  + data.cost[solut.s[i]][solut.s[i_next+1]];

        cost_new = solut.seq_get_C(0, i_prev)                                                    +           //       1st subseq
        solut.seq_get_W(i, i_next)               * (cost_concat_1) + data.cost[solut.s[i_next]][solut.s[i]]  +           // concat 2nd subseq
        solut.seq_get_W(i_next+1, data.dimen)    * (cost_concat_2) + solut.seq_get_C(i_next+1, data.dimen);   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data.dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut.seq_get_T(0, i_prev)       + data.cost[solut.s[i_prev]][solut.s[j]];
            cost_concat_2 = cost_concat_1                           + data.cost[solut.s[j]][solut.s[i_next]];
            cost_concat_3 = cost_concat_2 + solut.seq_get_T(i_next, j_prev)  + data.cost[solut.s[j_prev]][solut.s[i]];
            cost_concat_4 = cost_concat_3                           + data.cost[solut.s[i]][solut.s[j_next]];

            cost_new = solut.seq_get_C(0, i_prev)                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut.seq_get_W(i_next, j_prev)      * cost_concat_2 + solut.seq_get_C(i_next, j_prev) +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut.seq_get_W(j_next, data.dimen) * cost_concat_4 + solut.seq_get_C(j_next, data.dimen);   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        swap(solut.s, I, J);
        update_subseq_info_matrix(solut, data);

        assert(std::abs(cost_best-solut.cost) < DBL_EPSILON);
        assert(feasible(solut.s, data.dimen));
        return true;
    }

    return false;
}

bool search_two_opt(tSolution & solut, const tData & data) {
    Real cost_new,
        cost_concat_1, cost_concat_2;
    Real cost_best = DBL_MAX;
    Real rev_seq_cost;
    int i, j, i_prev, j_next;
    int I;
    int J;

    for (i = 1; i < data.dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut.seq_get_T(i, i+1);
        for (j = i + 2; j < data.dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += data.cost[solut.s[j-1]][solut.s[j]] * (solut.seq_get_W(i, j)-1.0);

          cost_concat_1 =                 solut.seq_get_T(0, i_prev)   + data.cost[solut.s[j]][solut.s[i_prev]];
          cost_concat_2 = cost_concat_1 + solut.seq_get_T(i, j)        + data.cost[solut.s[j_next]][solut.s[i]];

          cost_new = solut.seq_get_C(0, i_prev)                                                        +   //  1st subseq
              solut.seq_get_W(i, j)                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut.seq_get_W(j_next, data.dimen) * cost_concat_2 + solut.seq_get_C(j_next, data.dimen);      // concat 3rd subseq


            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        reverse(solut.s, I, J);
        update_subseq_info_matrix(solut, data);
        assert(std::abs(cost_best-solut.cost) < DBL_EPSILON);
        assert(feasible(solut.s, data.dimen));
        return true;
    }

    return false;
}

bool search_reinsertion(tSolution & solut, const tData & data, int opt) {
    Real cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    Real cost_best = DBL_MAX;
    int i, j, k, k_next, i_prev, j_next;
    int I;
    int J;
    int POS;

    for (i = 1, j = opt +i-1; i < data.dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 solut.seq_get_T(0, k)            + data.cost[solut.s[k]][solut.s[i]];
          cost_concat_2 = cost_concat_1 + solut.seq_get_T(i, j)            + data.cost[solut.s[j]][solut.s[k_next]];
          cost_concat_3 = cost_concat_2 + solut.seq_get_T(k_next, i_prev)  + data.cost[solut.s[i_prev]][solut.s[j_next]];

          cost_new = solut.seq_get_C(0, k)                                                                   +   //       1st subseq
              solut.seq_get_W(i, j)               * cost_concat_1 + solut.seq_get_C(i, j)                  +   //  concat 2nd subseq (reinserted seq)
              solut.seq_get_W(k_next, i_prev)     * cost_concat_2 + solut.seq_get_C(k_next, i_prev)        +   //  concat 3rd subseq
              solut.seq_get_W(j_next, data.dimen) * cost_concat_3 + solut.seq_get_C(j_next, data.dimen);       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data.dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut.seq_get_T(0, i_prev)  + data.cost[solut.s[i_prev]][solut.s[j_next]];
          cost_concat_2 = cost_concat_1 + solut.seq_get_T(j_next, k)  + data.cost[solut.s[k]][solut.s[i]];
          cost_concat_3 = cost_concat_2 + solut.seq_get_T(i, j)       + data.cost[solut.s[j]][solut.s[k_next]];

          cost_new = solut.seq_get_C(0, i_prev)                                                                  +   //       1st subseq
                  solut.seq_get_W(j_next, k)          * cost_concat_1 + solut.seq_get_C(j_next, k)             +   // concat 2nd subseq
                  solut.seq_get_W(i, j)               * cost_concat_2 + solut.seq_get_C(i, j)                  +   // concat 3rd subseq (reinserted seq)
                  solut.seq_get_W(k_next, data.dimen) * cost_concat_3 + solut.seq_get_C(k_next, data.dimen);       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        reinsert(solut.s, I, J, POS+1);
        update_subseq_info_matrix(solut, data);

        assert(std::abs(cost_best-solut.cost) < DBL_EPSILON);
        assert(feasible(solut.s, data.dimen));
        return true;
    }

    return false;
}


void RVND(tSolution & solut, tData & data) {

    int neighbd_list[] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    int nl_size = 5;
    int index;
    int neighbd;
    bool improve_flag;

    while (nl_size > 0) {
        index = data.rnd[data.rnd_index++];
        neighbd = neighbd_list[index];

        assert(index < nl_size);

        improve_flag = false;

        switch(neighbd){
            case REINSERTION:
                improve_flag = search_reinsertion(solut, data, REINSERTION);
                break;
            case OR_OPT_2:
                improve_flag = search_reinsertion(solut, data, OR_OPT_2);
                break;
            case OR_OPT_3:
                improve_flag = search_reinsertion(solut, data, OR_OPT_3);
                break;
            case SWAP:
                improve_flag = search_swap(solut, data);
                break;
            case TWO_OPT:
                improve_flag = search_two_opt(solut, data);
                break;
        }

        assert(feasible(solut.s, data.dimen));

        if (improve_flag) {
            neighbd_list[0] = SWAP;
            neighbd_list[1] = TWO_OPT;
            neighbd_list[2] = REINSERTION;
            neighbd_list[3] = OR_OPT_2;
            neighbd_list[4] = OR_OPT_3;
            nl_size = 5;
        } else {
            std::memmove(neighbd_list + index, neighbd_list + index+1, sizeof(int)*(nl_size-index-1));
            nl_size--;
        }


    }

}

void perturb(tSolution & solut_crnt, const tSolution & solut_partial, tData & data) {
    std::vector<int> s(solut_partial.s.begin(), solut_partial.s.begin()+(data.dimen+1));

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = floor((data.dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = data.rnd[data.rnd_index++];
        A_end = A_start + data.rnd[data.rnd_index++];

        B_start = data.rnd[data.rnd_index++];
        B_end = B_start + data.rnd[data.rnd_index++];
    }

    if (A_start < B_start) {
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    } else {
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }


    std::copy(s.begin(), s.begin()+(data.dimen+1), solut_crnt.s.begin());
}


void GILS_RVND(int Imax, int Iils, tData & data, int verbose) {

    tSolution solut_partial(data);
    tSolution solut_crnt(data);
    tSolution solut_best(data);

#ifdef MATRIX
    printf("MATRIX\n");
#elif defined(FLAT)
    printf("FLAT\n");
#endif
    printf("struct tSolution's size: %zu bytes\n", solut_crnt.size);
    printf("struct tSolution's size: %.4f MB\n", solut_crnt.MBsize);


    for(int i = 0; i < Imax; ++i){
        int aux = data.rnd[data.rnd_index++];

        Real alpha = R_table(aux);

        if (verbose) {
            printf("[+] Search %d\n", i+1);
            printf("\t[+] Constructing..\n");
        }


        construct(solut_crnt.s, alpha, data);
        assert(feasible(solut_crnt.s, data.dimen));
        if (verbose) {
            print_s(solut_crnt.s, data.dimen+1);
        }
        update_subseq_info_matrix(solut_crnt, data);

        solut_crnt.cpy_to(solut_partial, data);
        if (verbose) {
            printf("\t[+] Looking for the best Neighbor..\n");
            printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);
        }

        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(solut_crnt, data);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                solut_crnt.cpy_to(solut_partial, data);
                iterILS = 0;
            }

            perturb(solut_crnt, solut_partial, data);
            update_subseq_info_matrix(solut_crnt, data);

            assert(feasible(solut_crnt.s, data.dimen));
            iterILS++;
        }

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            solut_partial.cpy_to(solut_best, data);
        }

        if (verbose) {
            printf("\tCurrent best cost: %.2lf\n", solut_best.cost);

            printf("SOLUCAO: ");
            for(int i = 0; i < data.dimen+1; i++){
                printf("%d ", solut_best.s[i]);
            }
            printf("\n");
        }

    }
    printf("COST: %.2lf\n", solut_best.cost);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;
    int verbose = 0;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0) {
            verbose = 1;
        }
    }

    tData data;

    std::vector<int> rnd;

    data.dimen = loadData(data.cost, rnd);
    data.rnd = rnd;
    data.rnd_index = 0;

    Iils = data.dimen < 100 ? data.dimen : 100;
    std::clock_t start = clock();
    GILS_RVND(Imax, Iils, data, verbose);

    double res = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("TIME: %.6lf\n", res);

    return 0;
}
