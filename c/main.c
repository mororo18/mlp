#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "Data.h"

#define TRUE 1
#define FALSE 0

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

double ** cost;
int dimen;

typedef struct tRnd {
    int * rnd;
    int rnd_index;
} tRnd;

typedef struct tSubseq {
    double C;
    double T;
    double W;
} tSubseq;


double R_table(int i){
    static const double table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

void print_s(int * s, int sz) {

    for (int i = 0; i < sz; i++)
        printf("%d ", s[i]+1);
    printf("\n");
}


int partition(int *arr, int left, int right, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (cost[r][arr[j]] < cost[r][pivot]) {
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

void quicksort(int *arr, int left, int right, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, r);
        quicksort(arr, left, pivot - 1, r);
        quicksort(arr, pivot + 1, right, r);
    }
}

void sort(int *arr, int len, int r) {
    quicksort(arr, 0, len - 1, r);
}

void shift(int * vec, int from, int to, int sz) {
    if (from < to) {
        int dist = to - from;

        for (int i = from+sz-1, j = to+sz-1; i >= from; i--,j--) {
            vec[j] = vec[i];
        }
    } else {
        for (int i = from, j = to; i < from+sz; i++, j++) {
            vec[j] = vec[i];
        }
    }

}

void cpy(int * from, int * to, int sz) {
    for (int i = 0; i < sz; i++) {
        to[i] = from[i];
    }
}

void construct(int * ret, const double alpha, tRnd * rnd){

    int s[dimen+1];

    memset(s, 0, (dimen+1)*sizeof(int));

    int cL[dimen-1];
    int cL_size = dimen-1;
    
    for (int i = 1; i < dimen; ++i) {
        cL[i-1] = i;
    }

    int r = 0;
    for (int j = 1; j < dimen; ++j) {
        sort(cL, cL_size, r);

        int index = rnd->rnd[rnd->rnd_index++];
        int c = cL[index];
        s[j] = c;
        r = c;
        shift(cL, index+1, index, cL_size-index);
        cL_size--;
    }

    cpy(s, ret, (dimen+1));
}	

void swap(int * vec, int i, int j){
    int tmp = vec[i]; vec[i] = vec[j]; vec[j] = tmp;
}

void reverse(int * vec, int i, int j){
    int m = (i+j)/2;
    for (int first = i, last = j;
            first <= m;
            first++, last--) {

        swap(vec, first, last);
    }
}

void reinsert(int * vec, int i, int j, int pos){
    int seq[j-i+1];
    cpy(vec + i, seq, (j-i+1));

    if(pos < i){
        int sz = i-pos;
        shift(vec, pos, j+1-sz, sz);
        cpy(seq, vec + pos, j-i+1);
    }else{
        int sz = pos-j-1;
        shift(vec, j+1, i, sz);
        cpy(seq, vec+i+sz, j-i+1);
    }
}

double update_subseq_info_matrix(int * s, tSubseq seq[][dimen+1]){
    int i, j, j_prev, k;
    int from = 0;
    char t;
    for (i = 0; i < dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        seq[i][i].T = 0.0;
        seq[i][i].C = 0.0;
        seq[i][i].W = (double) !(i == 0);

        for (j = i+1; j < dimen+1; j++) {
            j_prev = j-1;
            
            seq[i][j].T = cost[s[j_prev]][s[j]] + seq[i][j_prev].T;
            seq[i][j].C = seq[i][j].T + seq[i][j_prev].C;
            seq[i][j].W = j + k;

        }
        from += t;
    }

    return seq[0][dimen].C;
}

double update_subseq_info_matrix_b(int * s, tSubseq seq[][dimen+1], int index){
    int i, j, j_prev, k;
    int from = index;
    char t;
    for (i = 0; i < dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        seq[i][i].T = 0.0;
        seq[i][i].C = 0.0;
        seq[i][i].W = (double) !(i == 0);

        for (j = i+1; j < dimen+1; j++) {
            j_prev = j-1;
            
            seq[i][j].T = cost[s[j_prev]][s[j]] + seq[i][j_prev].T;
            seq[i][j].C = seq[i][j].T + seq[i][j_prev].C;
            seq[i][j].W = j + k;

        }
        from += t;
    }

    return seq[0][dimen].C;
}

char search_swap(int * s, tSubseq seq[][dimen+1]) {
    double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    double cost_best = DBL_MAX;
    int i, j, j_prev, j_next, i_prev, i_next;
    int I;
    int J;

    for (i = 1; i < dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 seq[0][i_prev].T + cost[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[i][i_next].T  + cost[s[i]][s[i_next+1]];

        cost_new = seq[0][i_prev].C                                                    +           //       1st subseq
        seq[i][i_next].W               * (cost_concat_1) + cost[s[i_next]][s[i]]  +           // concat 2nd subseq
        seq[i_next+1][dimen].W   * (cost_concat_2) + seq[i_next+1][dimen].C;   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 seq[0][i_prev].T       + cost[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1                           + cost[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev].T  + cost[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3                           + cost[s[i]][s[j_next]];

            cost_new = seq[0][i_prev].C                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            seq[i_next][j_prev].W      * cost_concat_2 + seq[i_next][j_prev].C +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            seq[j_next][dimen].W * cost_concat_4 + seq[j_next][dimen].C;   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < seq[0][dimen].C - DBL_EPSILON) {
        swap(s, I, J);
        update_subseq_info_matrix_b(s, seq, I);
        return TRUE;
    }

    return FALSE;
}

char search_two_opt(int * s, tSubseq seq[][dimen+1]) {
    double cost_new, 
        cost_concat_1, cost_concat_2;
    double cost_best = DBL_MAX;// cost_l1, cost_l2;
    double rev_seq_cost;
    int i, j, i_prev, j_next;
    int I;
    int J;

    for (i = 1; i < dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = seq[i][i+1].T;
        for (j = i + 2; j < dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += cost[s[j-1]][s[j]] * (seq[i][j].W-1.0);

          cost_concat_1 =                 seq[0][i_prev].T   + cost[s[j]][s[i_prev]];
          cost_concat_2 = cost_concat_1 + seq[i][j].T        + cost[s[j_next]][s[i]];

          cost_new = seq[0][i_prev].C                                                        +   //  1st subseq
              seq[i][j].W                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              seq[j_next][dimen].W * cost_concat_2 + seq[j_next][dimen].C;      // concat 3rd subseq

            
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < seq[0][dimen].C - DBL_EPSILON) {
        reverse(s, I, J);
        update_subseq_info_matrix(s, seq);
        return TRUE;
    }

    return FALSE;
}

char search_reinsertion(int * s, tSubseq seq[][dimen+1], const int opt) {
    double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    int i, j, k, k_next, i_prev, j_next;
    int I;
    int J;
    int POS;

    for (i = 1, j = opt +i-1; i < dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 seq[0][k].T            + cost[s[k]][s[i]];
          cost_concat_2 = cost_concat_1 + seq[i][j].T            + cost[s[j]][s[k_next]];
          cost_concat_3 = cost_concat_2 + seq[k_next][i_prev].T  + cost[s[i_prev]][s[j_next]];

          cost_new = seq[0][k].C                                                                   +   //       1st subseq
              seq[i][j].W               * cost_concat_1 + seq[i][j].C                  +   //  concat 2nd subseq (reinserted seq)
              seq[k_next][i_prev].W     * cost_concat_2 + seq[k_next][i_prev].C        +   //  concat 3rd subseq
              seq[j_next][dimen].W * cost_concat_3 + seq[j_next][dimen].C;       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 seq[0][i_prev].T  + cost[s[i_prev]][s[j_next]];
          cost_concat_2 = cost_concat_1 + seq[j_next][k].T  + cost[s[k]][s[i]];
          cost_concat_3 = cost_concat_2 + seq[i][j].T       + cost[s[j]][s[k_next]];

          cost_new = seq[0][i_prev].C                                                                  +   //       1st subseq
                  seq[j_next][k].W          * cost_concat_1 + seq[j_next][k].C             +   // concat 2nd subseq
                  seq[i][j].W               * cost_concat_2 + seq[i][j].C                  +   // concat 3rd subseq (reinserted seq)
                  seq[k_next][dimen].W * cost_concat_3 + seq[k_next][dimen].C;       // concat 4th subseq
          
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < seq[0][dimen].C - DBL_EPSILON) {
        reinsert(s, I, J, POS+1);
        update_subseq_info_matrix_b(s, seq, I < POS+1 ? I : POS+1);
        return TRUE;
    }

    return FALSE;
}


void RVND(int * solut, tSubseq seq[][dimen+1], tRnd * rnd) {

    int neighbd_list[] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    int nl_size = 5;
    int index;
    int neighbd;
    char improve_flag;

    while (nl_size > 0) {

        index = rnd->rnd[rnd->rnd_index++];
        neighbd = neighbd_list[index];

        improve_flag = FALSE;

        switch(neighbd){
            case REINSERTION:
                improve_flag = search_reinsertion(solut, seq, REINSERTION);
                break;				
            case OR_OPT_2:
                improve_flag = search_reinsertion(solut, seq, OR_OPT_2);
                break;				
            case OR_OPT_3:
                improve_flag = search_reinsertion(solut, seq, OR_OPT_3);
                break;				
            case SWAP:
                improve_flag = search_swap(solut, seq);
                break;
            case TWO_OPT:
                improve_flag = search_two_opt(solut, seq);
                break;				
        }

        if (improve_flag) {
            neighbd_list[0] = SWAP;
            neighbd_list[1] = TWO_OPT;
            neighbd_list[2] = REINSERTION;
            neighbd_list[3] = OR_OPT_2;
            neighbd_list[4] = OR_OPT_3;
            nl_size = 5;
        } else {
            shift(neighbd_list, index+1, index, (nl_size-index-1));
            nl_size--;
        }

    }

}

void perturb(int * s_crnt, int * s_partial, tRnd * rnd) {
    int s[dimen+1];
    cpy(s_partial, s, (dimen+1));

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        A_start = rnd->rnd[rnd->rnd_index++];
        A_end = A_start + rnd->rnd[rnd->rnd_index++];

        B_start = rnd->rnd[rnd->rnd_index++];
        B_end = B_start + rnd->rnd[rnd->rnd_index++];
    }
    
    if (A_start < B_start) {
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    } else {
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    cpy(s, s_crnt, (dimen+1));
}


void GILS_RVND(int Imax, int Iils, tRnd * rnd) {

    int s_partial[dimen+1];
    int s_crnt[dimen+1];
    int s_best[dimen+1];

    double cost_best = DBL_MAX;
    double cost_crnt;
    double cost_partial;

    tSubseq seq[dimen+1][dimen+1];

    memset(seq, 0, sizeof(tSubseq)*(dimen+1)*(dimen+1));
    memset(s_partial, 0, sizeof(int)*(dimen+1));
    memset(s_crnt, 0, sizeof(int)*(dimen+1));
    memset(s_best, 0, sizeof(int)*(dimen+1));

    for(int i = 0; i < Imax; ++i){
        int aux = rnd->rnd[rnd->rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	

        construct(s_crnt, alpha, rnd);
        print_s(s_crnt, dimen+1);
        cost_crnt = update_subseq_info_matrix(s_crnt, seq);

        cpy(s_crnt, s_partial, (dimen+1));
        cost_partial = cost_crnt;

        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", cost_partial);	

        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(s_crnt, seq, rnd);
            cost_crnt = seq[0][dimen].C;
            if(cost_crnt < cost_partial - DBL_EPSILON){
                cpy(s_crnt, s_partial, (dimen+1));
                cost_partial = cost_crnt;
                iterILS = 0;
            }

            perturb(s_crnt, s_partial, rnd);
            update_subseq_info_matrix(s_crnt, seq);
            iterILS++;
        }

        if (cost_partial < cost_best - DBL_EPSILON) {
            cpy(s_partial, s_best, (dimen+1));
            cost_best = cost_partial;
        }

        printf("\tCurrent best cost: %.2lf\n", cost_best);
        printf("SOLUCAO: ");
        for(int i = 0; i < dimen+1; i++){
            printf("%d ", s_best[i]);
        }
        printf("\n");

    }
    printf("COST: %.2lf\n", cost_best);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tRnd rnd = {0};

    int * rnd_adr;

    dimen = loadData(&cost, &rnd_adr);
    rnd.rnd = rnd_adr;
    rnd.rnd_index = 0;

    srand(clock());

    Iils = dimen < 100 ? dimen : 100;
    time_t start = clock();
    GILS_RVND(Imax, Iils, &rnd);

    double res = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("TIME: %.6lf\n", res);

    return 0;
}

