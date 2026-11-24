#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <stdalign.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "readData.h"
#include "Data.h"

#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

#define FLAT

typedef unsigned uint;

typedef struct tData {
    double ** cost;
    int dimen;
    int * rnd;
    int rnd_size;
    uint rnd_index;
} tData;

typedef struct tSeqInfo {
    double T;
    double C;
    double W;
} tSeqInfo;

enum Neighborhood {
    REINSERTION = 1,
    OR_OPT_2,
    OR_OPT_3,
    SWAP,
    TWO_OPT
};

typedef struct tSolution {
    int * s;
    int s_size;
#ifdef FLAT
    tSeqInfo * seq;
#elif defined(MATRIX)
    tSeqInfo ** seq;
#endif
    double cost;
} tSolution;

static inline unsigned to_1D(int x, int y, int size) {
    int a = x*(size+1) + y;

    return a;
}

tSolution Solution_init(tData data) {
    tSolution solut;
    solut.s = (int *) calloc((size_t)(data.dimen+1), sizeof(int));
    solut.s_size = data.dimen+1;

#ifdef FLAT
    solut.seq = (tSeqInfo *) malloc((size_t)(data.dimen+1) * (data.dimen+1) * sizeof(tSeqInfo));
#elif defined(MATRIX)

    solut.seq = (tSeqInfo **) malloc((size_t)(data.dimen+1) * sizeof(tSeqInfo *));
    for (int i = 0; i < data.dimen+1; i++) {
        solut.seq[i] = (tSeqInfo *) malloc((size_t)(data.dimen+1) * sizeof(tSeqInfo));
        memset(solut.seq[i], 0, (data.dimen+1)*sizeof(tSeqInfo));
    }

#endif
    solut.cost = DBL_MAX;

    return solut;
}

void Solution_cpy(tSolution * src, tSolution * tgt, const tData * data) {

    memcpy(tgt->s, src->s, (size_t) src->s_size * sizeof(int));
    tgt->s_size = src->s_size;
    tgt->cost = src->cost;

}

double R_table(int i){
    static const double table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

void print_s(int * s, int size) {

    for (int i = 0; i < size; i++)
        printf("%d ", s[i]+1);
    printf("\n");
}

static inline void swap_int(int * a, int * b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

int partition(int * arr, int left, int right, const tData * data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (data->cost[r][arr[j]] < data->cost[r][pivot]) {
            i++;
            swap_int(&arr[i], &arr[j]);
        }
    }
    swap_int(&arr[i + 1], &arr[right]);
    return i + 1;
}

void quicksort(int * arr, int left, int right, const tData * data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

void sort(int * arr, int size, int r, const tData * data) {
    quicksort(arr, 0, size - 1, data, r);
}

int construct(const double alpha, tData * data, int * s){

    int s_size = 0;
    s[s_size++] = 0;

    int cL_size = data->dimen-1;
    int * cL = (int *) malloc((size_t) (data->dimen-1) * sizeof(int));
    for(int i = 1; i < data->dimen; ++i){
        cL[i-1] = i;
    }

    int r = 0;
    while (cL_size > 0) {
        sort(cL, cL_size, r, data);

        int index = data->rnd[data->rnd_index++];

        int c = cL[index];
        s[s_size++] = c;
        r = c;
        memmove(cL + index, cL + index + 1, (size_t)(cL_size - index - 1) * sizeof(int));
        cL_size--;
    }

    s[s_size++] = 0;

    free(cL);

    return s_size;
}

void swap(int * vec, int i, int j){
    swap_int(&vec[i], &vec[j]);
}

void reverse(int * vec, int i, int j){
    while (i < j) {
        swap_int(&vec[i], &vec[j]);
        i++;
        j--;
    }
}

void reinsert(int * vec, int i, int j, int pos){
    int sz = j - i + 1;
    int * seq = (int *) malloc((size_t) sz * sizeof(int));
    memcpy(seq, vec + i, (size_t) sz * sizeof(int));

    if(pos < i){
        memmove(vec + pos + sz, vec + pos, (size_t)(i - pos) * sizeof(int));
        memcpy(vec + pos, seq, (size_t) sz * sizeof(int));
    }else{
        memmove(vec + i, vec + j + 1, (size_t)(pos - j - 1) * sizeof(int));
        memcpy(vec + pos - sz, seq, (size_t) sz * sizeof(int));
    }

    free(seq);
}

#ifdef FLAT

void update_subseq_info_matrix(tSolution * solut, const tData * data, int index){
    alignas(INT_SZ) int i, j, j_prev, k;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for (i = 0; i < data->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut->seq[to_1D(i, i, data->dimen)].T = 0.0;
        solut->seq[to_1D(i, i, data->dimen)].C = 0.0;
        solut->seq[to_1D(i, i, data->dimen)].W = (double) !(i == 0);

        for (j = i+1; j < data->dimen+1; j++) {
            j_prev = j-1;

            solut->seq[to_1D(i, j, data->dimen)].T = data->cost[solut->s[j_prev]][solut->s[j]] + solut->seq[to_1D(i, j_prev, data->dimen)].T;
            solut->seq[to_1D(i, j, data->dimen)].C = solut->seq[to_1D(i, j, data->dimen)].T + solut->seq[to_1D(i, j_prev, data->dimen)].C;
            solut->seq[to_1D(i, j, data->dimen)].W = j + k;
        }
        from += t;
    }

    solut->cost = solut->seq[to_1D(0, data->dimen, data->dimen)].C;
}

bool search_swap(tSolution * solut, const tData * data) {
    alignas(DBL_SZ) double cost_new,
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    alignas(DBL_SZ) double cost_best = DBL_MAX;
    alignas(INT_SZ) int i, j, j_prev, j_next, i_prev, i_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut->seq[to_1D(0, i_prev, data->dimen)].T + data->cost[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + solut->seq[to_1D(i, i_next, data->dimen)].T  + data->cost[solut->s[i]][solut->s[i_next+1]];

        cost_new = solut->seq[to_1D(0, i_prev, data->dimen)].C                                                    +           //       1st subseq
        solut->seq[to_1D(i, i_next, data->dimen)].W               * (cost_concat_1) + data->cost[solut->s[i_next]][solut->s[i]]  +           // concat 2nd subseq
        solut->seq[to_1D(i_next+1, data->dimen, data->dimen)].W   * (cost_concat_2) + solut->seq[to_1D(i_next+1, data->dimen, data->dimen)].C;   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data->dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut->seq[to_1D(0, i_prev, data->dimen)].T       + data->cost[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1                           + data->cost[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + solut->seq[to_1D(i_next, j_prev, data->dimen)].T  + data->cost[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3                           + data->cost[solut->s[i]][solut->s[j_next]];

            cost_new = solut->seq[to_1D(0, i_prev, data->dimen)].C                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut->seq[to_1D(i_next, j_prev, data->dimen)].W      * cost_concat_2 + solut->seq[to_1D(i_next, j_prev, data->dimen)].C +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut->seq[to_1D(j_next, data->dimen, data->dimen)].W * cost_concat_4 + solut->seq[to_1D(j_next, data->dimen, data->dimen)].C;   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        swap(solut->s, I, J);
        update_subseq_info_matrix(solut, data, I);
        return true;
    }

    return false;
}

bool search_two_opt(tSolution * solut, const tData * data) {
    alignas(DBL_SZ) double cost_new,
        cost_concat_1, cost_concat_2;
    alignas(DBL_SZ) double cost_best = DBL_MAX; // cost_l1, cost_l2;
    alignas(DBL_SZ) double rev_seq_cost;
    alignas(INT_SZ) int i, j, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut->seq[to_1D(i, i+1, data->dimen)].T;
        for (j = i + 2; j < data->dimen; ++j) {
            j_next = j + 1;

            rev_seq_cost += data->cost[solut->s[j-1]][solut->s[j]] * (solut->seq[to_1D(i, j, data->dimen)].W-1.0);

            cost_concat_1 =                 solut->seq[to_1D(0, i_prev, data->dimen)].T   + data->cost[solut->s[j]][solut->s[i_prev]];
            cost_concat_2 = cost_concat_1 + solut->seq[to_1D(i, j, data->dimen)].T        + data->cost[solut->s[j_next]][solut->s[i]];

            cost_new = solut->seq[to_1D(0, i_prev, data->dimen)].C                                                        +   //  1st subseq
                    solut->seq[to_1D(i, j, data->dimen)].W                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
                    solut->seq[to_1D(j_next, data->dimen, data->dimen)].W * cost_concat_2 + solut->seq[to_1D(j_next, data->dimen, data->dimen)].C;      // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reverse(solut->s, I, J);
        update_subseq_info_matrix(solut, data, 0);
        return true;
    }

    return false;
}

bool search_reinsertion(tSolution * solut, const tData * data, const int opt) {
    alignas(DBL_SZ) double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    alignas(DBL_SZ) double cost_best = DBL_MAX; //, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, k_next, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;
    alignas(INT_SZ) int POS;

    for (i = 1, j = opt +i-1; i < data->dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

            cost_concat_1 =                 solut->seq[to_1D(0, k, data->dimen)].T            + data->cost[solut->s[k]][solut->s[i]];
            cost_concat_2 = cost_concat_1 + solut->seq[to_1D(i, j, data->dimen)].T            + data->cost[solut->s[j]][solut->s[k_next]];
            cost_concat_3 = cost_concat_2 + solut->seq[to_1D(k_next, i_prev, data->dimen)].T  + data->cost[solut->s[i_prev]][solut->s[j_next]];

            cost_new = solut->seq[to_1D(0, k, data->dimen)].C                                                                   +   //       1st subseq
                solut->seq[to_1D(i, j, data->dimen)].W               * cost_concat_1 + solut->seq[to_1D(i, j, data->dimen)].C                  +   //  concat 2nd subseq (reinserted seq)
                solut->seq[to_1D(k_next, i_prev, data->dimen)].W     * cost_concat_2 + solut->seq[to_1D(k_next, i_prev, data->dimen)].C        +   //  concat 3rd subseq
                solut->seq[to_1D(j_next, data->dimen, data->dimen)].W * cost_concat_3 + solut->seq[to_1D(j_next, data->dimen, data->dimen)].C;       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data->dimen; ++k) {
            k_next = k + 1;

            cost_concat_1 =                 solut->seq[to_1D(0, i_prev, data->dimen)].T  + data->cost[solut->s[i_prev]][solut->s[j_next]];
            cost_concat_2 = cost_concat_1 + solut->seq[to_1D(j_next, k, data->dimen)].T  + data->cost[solut->s[k]][solut->s[i]];
            cost_concat_3 = cost_concat_2 + solut->seq[to_1D(i, j, data->dimen)].T       + data->cost[solut->s[j]][solut->s[k_next]];

            cost_new = solut->seq[to_1D(0, i_prev, data->dimen)].C                                                                  +   //       1st subseq
                solut->seq[to_1D(j_next, k, data->dimen)].W          * cost_concat_1 + solut->seq[to_1D(j_next, k, data->dimen)].C             +   // concat 2nd subseq
                solut->seq[to_1D(i, j, data->dimen)].W               * cost_concat_2 + solut->seq[to_1D(i, j, data->dimen)].C                  +   // concat 3rd subseq (reinserted seq)
                solut->seq[to_1D(k_next, data->dimen, data->dimen)].W * cost_concat_3 + solut->seq[to_1D(k_next, data->dimen, data->dimen)].C;       // concat 4th subseq


            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reinsert(solut->s, I, J, POS+1);
        update_subseq_info_matrix(solut, data, I < POS+1 ? I : POS+1);
        return true;
    }

    return false;
}

#elif defined(MATRIX)

void update_subseq_info_matrix(tSolution * solut, const tData * data, int index){
    alignas(INT_SZ) int i, j, j_prev, k;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for (i = 0; i < data->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut->seq[i][i].T = 0.0;
        solut->seq[i][i].C = 0.0;
        solut->seq[i][i].W = (double) !(i == 0);

        for (j = i+1; j < data->dimen+1; j++) {
            j_prev = j-1;

            solut->seq[i][j].T = data->cost[solut->s[j_prev]][solut->s[j]] + solut->seq[i][j_prev].T;
            solut->seq[i][j].C = solut->seq[i][j].T + solut->seq[i][j_prev].C;
            solut->seq[i][j].W = j + k;

        }
        from += t;
    }

    solut->cost = solut->seq[0][data->dimen].C;
}

bool search_swap(tSolution * solut, const tData * data) {
    alignas(DBL_SZ) double cost_new,
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    alignas(DBL_SZ) double cost_best = DBL_MAX;
    alignas(INT_SZ) int i, j, j_prev, j_next, i_prev, i_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut->seq[0][i_prev].T + data->cost[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + solut->seq[i][i_next].T  + data->cost[solut->s[i]][solut->s[i_next+1]];

        cost_new = solut->seq[0][i_prev].C                                                    +           //       1st subseq
        solut->seq[i][i_next].W               * (cost_concat_1) + data->cost[solut->s[i_next]][solut->s[i]]  +           // concat 2nd subseq
        solut->seq[i_next+1][data->dimen].W   * (cost_concat_2) + solut->seq[i_next+1][data->dimen].C;   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data->dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut->seq[0][i_prev].T       + data->cost[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1                           + data->cost[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + solut->seq[i_next][j_prev].T  + data->cost[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3                           + data->cost[solut->s[i]][solut->s[j_next]];

            cost_new = solut->seq[0][i_prev].C                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut->seq[i_next][j_prev].W      * cost_concat_2 + solut->seq[i_next][j_prev].C +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut->seq[j_next][data->dimen].W * cost_concat_4 + solut->seq[j_next][data->dimen].C;   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        swap(solut->s, I, J);
        update_subseq_info_matrix(solut, data, I);
        return true;
    }

    return false;
}

bool search_two_opt(tSolution * solut, const tData * data) {
    alignas(DBL_SZ) double cost_new,
        cost_concat_1, cost_concat_2;
    alignas(DBL_SZ) double cost_best = DBL_MAX; // cost_l1, cost_l2;
    alignas(DBL_SZ) double rev_seq_cost;
    alignas(INT_SZ) int i, j, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut->seq[i][i+1].T;
        for (j = i + 2; j < data->dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += data->cost[solut->s[j-1]][solut->s[j]] * (solut->seq[i][j].W-1.0);

          cost_concat_1 =                 solut->seq[0][i_prev].T   + data->cost[solut->s[j]][solut->s[i_prev]];
          cost_concat_2 = cost_concat_1 + solut->seq[i][j].T        + data->cost[solut->s[j_next]][solut->s[i]];

          cost_new = solut->seq[0][i_prev].C                                                        +   //  1st subseq
              solut->seq[i][j].W                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut->seq[j_next][data->dimen].W * cost_concat_2 + solut->seq[j_next][data->dimen].C;      // concat 3rd subseq


            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reverse(solut->s, I, J);
        update_subseq_info_matrix(solut, data, 0);
        return true;
    }

    return false;
}

bool search_reinsertion(tSolution * solut, const tData * data, const int opt) {
    alignas(DBL_SZ) double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    alignas(DBL_SZ) double cost_best = DBL_MAX; //, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, k_next, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;
    alignas(INT_SZ) int POS;

    for (i = 1, j = opt +i-1; i < data->dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 solut->seq[0][k].T            + data->cost[solut->s[k]][solut->s[i]];
          cost_concat_2 = cost_concat_1 + solut->seq[i][j].T            + data->cost[solut->s[j]][solut->s[k_next]];
          cost_concat_3 = cost_concat_2 + solut->seq[k_next][i_prev].T  + data->cost[solut->s[i_prev]][solut->s[j_next]];

          cost_new = solut->seq[0][k].C                                                                   +   //       1st subseq
              solut->seq[i][j].W               * cost_concat_1 + solut->seq[i][j].C                  +   //  concat 2nd subseq (reinserted seq)
              solut->seq[k_next][i_prev].W     * cost_concat_2 + solut->seq[k_next][i_prev].C        +   //  concat 3rd subseq
              solut->seq[j_next][data->dimen].W * cost_concat_3 + solut->seq[j_next][data->dimen].C;       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data->dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut->seq[0][i_prev].T  + data->cost[solut->s[i_prev]][solut->s[j_next]];
          cost_concat_2 = cost_concat_1 + solut->seq[j_next][k].T  + data->cost[solut->s[k]][solut->s[i]];
          cost_concat_3 = cost_concat_2 + solut->seq[i][j].T       + data->cost[solut->s[j]][solut->s[k_next]];

          cost_new = solut->seq[0][i_prev].C                                                                  +   //       1st subseq
                  solut->seq[j_next][k].W          * cost_concat_1 + solut->seq[j_next][k].C             +   // concat 2nd subseq
                  solut->seq[i][j].W               * cost_concat_2 + solut->seq[i][j].C                  +   // concat 3rd subseq (reinserted seq)
                  solut->seq[k_next][data->dimen].W * cost_concat_3 + solut->seq[k_next][data->dimen].C;       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reinsert(solut->s, I, J, POS+1);
        update_subseq_info_matrix(solut, data, I < POS+1 ? I : POS+1);
        return true;
    }

    return false;
}

#endif

void RVND(tSolution * solut, tData * data) {

    alignas(INT_SZ) int neighbd_list[5] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    alignas(INT_SZ) int neighbd_list_size = 5;
    alignas(INT_SZ) uint index;
    alignas(INT_SZ) int neighbd;
    bool improve_flag;

    while (neighbd_list_size > 0) {

        index = data->rnd[data->rnd_index++];
        neighbd = neighbd_list[index];
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
        if (improve_flag) {
            neighbd_list[0] = SWAP;
            neighbd_list[1] = TWO_OPT;
            neighbd_list[2] = REINSERTION;
            neighbd_list[3] = OR_OPT_2;
            neighbd_list[4] = OR_OPT_3;
            neighbd_list_size = 5;
        } else {
            memmove(neighbd_list + index, neighbd_list + index + 1, (size_t)(neighbd_list_size - index - 1) * sizeof(int));
            neighbd_list_size--;
        }

    }

}

void perturb(tSolution * solut, tData * data, int * s) {
    memcpy(s, solut->s, (size_t) solut->s_size * sizeof(int));
    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = (int) floor((data->dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;
    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = data->rnd[data->rnd_index++];
        A_end = A_start + data->rnd[data->rnd_index++];

        B_start = data->rnd[data->rnd_index++];
        B_end = B_start + data->rnd[data->rnd_index++];
    }

    if (A_start < B_start) {
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    } else {
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    (void) size_max;
    (void) size_min;
}


void GILS_RVND(int Imax, int Iils, tData * data, bool verbose) {

    tSolution solut_partial = Solution_init(*data);
    tSolution solut_crnt = Solution_init(*data);
    tSolution solut_best = Solution_init(*data);

    for(int i = 0; i < Imax; ++i){
        int aux = data->rnd[data->rnd_index++];

        double alpha = R_table(aux);

        if (verbose) {
            printf("[+] Search %d\n", i+1);
            printf("\t[+] Constructing..\n");
        }


        solut_crnt.s_size = construct(alpha, data, solut_crnt.s);
        update_subseq_info_matrix(&solut_crnt, data, 0);

        Solution_cpy(&solut_crnt, &solut_partial, data);
        if (verbose) {
            printf("\t[+] Looking for the best Neighbor..\n");
            printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);
        }

        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(&solut_crnt, data);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                Solution_cpy(&solut_crnt, &solut_partial, data);
                iterILS = 0;
            }

            perturb(&solut_partial, data, solut_crnt.s);
            solut_crnt.s_size = solut_partial.s_size;
            update_subseq_info_matrix(&solut_crnt, data, 0);
            iterILS++;
        }

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(&solut_partial, &solut_best, data);
        }

        if (verbose) {
            printf("\tCurrent best cost: %g\n", solut_best.cost);

            printf("SOLUCAO: ");
            for(int j = 0; j < solut_best.s_size; j++){
                printf("%d ", solut_best.s[j]);
            }
            printf("\n");
        }

    }
    printf("COST: %.2lf\n", solut_best.cost);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;
    bool verbose = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    tData data;
    memset(&data, 0, sizeof(data));


    int * rnd = NULL;
    int rnd_size = 0;

    data.dimen = loadData(&data.cost, &rnd, &rnd_size);
    data.rnd = rnd;
    data.rnd_size = rnd_size;
    data.rnd_index = 0;

    tSolution solut = Solution_init(data);
    (void) solut;

    srand(clock());

    Iils = data.dimen < 100 ? data.dimen : 100;
    struct timespec t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    size_t start = clock();
    GILS_RVND(Imax, Iils, &data, verbose);
    double cpu_time = (double)(clock() - start) / CLOCKS_PER_SEC ;
    clock_gettime(CLOCK_MONOTONIC, &t2);
    long duration = (t2.tv_sec - t1.tv_sec) * 1000L + (t2.tv_nsec - t1.tv_nsec) / 1000000L;

    double res = (double)duration / 10e2;
    printf("TIME: %g\n", cpu_time);
    printf("wall clock (s): %g\n", res);

    printf("Tamanho RND %d\n", rnd_size);

    exit(0);

    return 0;
}
