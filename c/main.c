#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "data.h"
#include "types.h"

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

Real R_table(int i){
    static const Real table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

void print_s(int * s, int sz) {

    for (int i = 0; i < sz; i++)
        printf("%d ", s[i]);
    printf("\n");
}

int partition(int *arr, int left, int right, tData *data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (data->cost[r][arr[j]] < data->cost[r][pivot]) {
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

void quicksort(int *arr, int left, int right, tData *data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

void sort(int *arr, int len, int r, tData *data) {
    quicksort(arr, 0, len - 1, data, r);
}

#ifndef NDEBUG
_Bool feasible(int * s, int sz) {
    assert(sz > 0);
    int * count = calloc(sz, sizeof(int)); 

    for (int i = 0; i < sz; i++) { 
        assert(s[i] >= 0 && s[i] < sz);
    }

    for (int i = 0; i < sz; i++) 
        count[s[i]] += 1;

    for (int i = 0; i < sz; i++)
        if (count[i] != 1) {
            free(count);
            return false;
        }

    free(count);
    return true;
}
#endif

void construct(int * ret, const Real alpha, tData * data){

    int * s = calloc(data->dimen+1, sizeof(int));
    int * cL = calloc(data->dimen-1, sizeof(int));
    int cL_size = data->dimen-1;
    
    for (int i = 1; i < data->dimen; ++i) {
        cL[i-1] = i;
    }

    int r = 0;
    for (int j = 1; j < data->dimen; ++j) {
        sort(cL, cL_size, r, data);

        int index = data->rnd[data->rnd_index++];
        int c = cL[index];
        s[j] = c;
        r = c;
        memmove(cL+index, cL+index+1, sizeof(int)*(cL_size-index));
        cL_size--;
    }

    memcpy(ret, s, sizeof(int)*(data->dimen+1));
    free(s);
    free(cL);
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
    memcpy(seq, vec + i, sizeof(int)*(j-i+1));

    if(pos < i){
        int sz = i-pos;
        memmove(vec + j+1-sz, vec + pos, sizeof(int)*sz);
        memcpy(vec + pos, seq, sizeof(int)*(j-i+1));
    }else{
        int sz = pos-j-1;
        memmove(vec + i, vec + j+1, sizeof(int)*sz);
        memcpy(vec + i+sz, seq, sizeof(int)*(j-i+1));
    }
}

void update_subseq_info_matrix(tSolution * solut, const tData * data){
    int i, j, j_prev, k;
    int from = 0;
    _Bool t;
    for (i = 0; i < data->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        seq_set_T(solut, i, i, 0.0);
        seq_set_C(solut, i, i, 0.0);
        seq_set_W(solut, i, i, (Real)!(i == 0));

        for (j = i+1; j < data->dimen+1; j++) {
            j_prev = j-1;
            
            Real T = data->cost[solut->s[j_prev]][solut->s[j]] + seq_get_T(solut,  i,  j_prev);
            seq_set_T(solut, i, j, T);

            Real C = seq_get_T(solut,  i,  j) + seq_get_C(solut,  i,  j_prev);
            seq_set_C(solut, i, j, C);

            Real W = j + k;
            seq_set_W(solut, i, j, W);

        }
        from += t;
    }

    solut->cost = seq_get_C(solut,  0,  data->dimen);
}

void update_subseq_info_matrix_b(tSolution * solut, const tData * data, int index){
    int i, j, j_prev, k;
    int from = index;
    _Bool t;
    for (i = 0; i < data->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        seq_set_T(solut, i, i, 0.0);
        seq_set_C(solut, i, i, 0.0);
        seq_set_W(solut, i, i, (Real)!(i == 0));

        for (j = i+1; j < data->dimen+1; j++) {
            j_prev = j-1;
            
            Real T = data->cost[solut->s[j_prev]][solut->s[j]] + seq_get_T(solut,  i,  j_prev);
            seq_set_T(solut, i, j, T);

            Real C = seq_get_T(solut,  i,  j) + seq_get_C(solut,  i,  j_prev);
            seq_set_C(solut, i, j, C);

            Real W = j + k;
            seq_set_W(solut, i, j, W);


        }
        from += t;
    }

    solut->cost = seq_get_C(solut,  0,  data->dimen);
}

_Bool search_swap(tSolution * solut, const tData * data) {
    Real cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    Real cost_best = DBL_MAX;
    int i, j, j_prev, j_next, i_prev, i_next;
    int I;
    int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 seq_get_T(solut,  0,  i_prev) + data->cost[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + seq_get_T(solut,  i,  i_next)  + data->cost[solut->s[i]][solut->s[i_next+1]];

        cost_new = seq_get_C(solut,  0,  i_prev)                                                    +           //       1st subseq
        seq_get_W(solut,  i,  i_next)               * (cost_concat_1) + data->cost[solut->s[i_next]][solut->s[i]]  +           // concat 2nd subseq
        seq_get_W(solut,  i_next+1,  data->dimen)   * (cost_concat_2) + seq_get_C(solut,  i_next+1,  data->dimen);   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data->dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 seq_get_T(solut,  0,  i_prev)       + data->cost[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1                           + data->cost[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + seq_get_T(solut,  i_next,  j_prev)  + data->cost[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3                           + data->cost[solut->s[i]][solut->s[j_next]];

            cost_new = seq_get_C(solut,  0,  i_prev)                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            seq_get_W(solut,  i_next,  j_prev)      * cost_concat_2 + seq_get_C(solut,  i_next,  j_prev) +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            seq_get_W(solut,  j_next,  data->dimen) * cost_concat_4 + seq_get_C(solut,  j_next,  data->dimen);   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        swap(solut->s, I, J);
        update_subseq_info_matrix_b(solut, data, I);

        assert(abs(cost_best-solut->cost) < DBL_EPSILON);
        assert(feasible(solut->s, data->dimen));
        return true;
    }

    return false;
}

_Bool search_two_opt(tSolution * solut, const tData * data) {
    Real cost_new, 
        cost_concat_1, cost_concat_2;
    Real cost_best = DBL_MAX;// cost_l1, cost_l2;
    Real rev_seq_cost;
    int i, j, i_prev, j_next;
    int I;
    int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = seq_get_T(solut,  i,  i+1);
        for (j = i + 2; j < data->dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += data->cost[solut->s[j-1]][solut->s[j]] * (seq_get_W(solut,  i,  j)-1.0);

          cost_concat_1 =                 seq_get_T(solut,  0,  i_prev)   + data->cost[solut->s[j]][solut->s[i_prev]];
          cost_concat_2 = cost_concat_1 + seq_get_T(solut,  i,  j)        + data->cost[solut->s[j_next]][solut->s[i]];

          cost_new = seq_get_C(solut,  0,  i_prev)                                                        +   //  1st subseq
              seq_get_W(solut,  i,  j)                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              seq_get_W(solut,  j_next,  data->dimen) * cost_concat_2 + seq_get_C(solut,  j_next,  data->dimen);      // concat 3rd subseq

            
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reverse(solut->s, I, J);
        update_subseq_info_matrix(solut, data);
        assert(abs(cost_best-solut->cost) < DBL_EPSILON);
        assert(feasible(solut->s, data->dimen));
        return true;
    }

    return false;
}

_Bool search_reinsertion(tSolution * solut, const tData * data, int opt) {
    Real cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    Real cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    int i, j, k, k_next, i_prev, j_next;
    int I;
    int J;
    int POS;

    for (i = 1, j = opt +i-1; i < data->dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 seq_get_T(solut,  0,  k)            + data->cost[solut->s[k]][solut->s[i]];
          cost_concat_2 = cost_concat_1 + seq_get_T(solut,  i,  j)            + data->cost[solut->s[j]][solut->s[k_next]];
          cost_concat_3 = cost_concat_2 + seq_get_T(solut,  k_next,  i_prev)  + data->cost[solut->s[i_prev]][solut->s[j_next]];

          cost_new = seq_get_C(solut,  0,  k)                                                                   +   //       1st subseq
              seq_get_W(solut,  i,  j)               * cost_concat_1 + seq_get_C(solut,  i,  j)                  +   //  concat 2nd subseq (reinserted seq)
              seq_get_W(solut,  k_next,  i_prev)     * cost_concat_2 + seq_get_C(solut,  k_next,  i_prev)        +   //  concat 3rd subseq
              seq_get_W(solut,  j_next,  data->dimen) * cost_concat_3 + seq_get_C(solut,  j_next,  data->dimen);       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data->dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 seq_get_T(solut,  0,  i_prev)  + data->cost[solut->s[i_prev]][solut->s[j_next]];
          cost_concat_2 = cost_concat_1 + seq_get_T(solut,  j_next,  k)  + data->cost[solut->s[k]][solut->s[i]];
          cost_concat_3 = cost_concat_2 + seq_get_T(solut,  i,  j)       + data->cost[solut->s[j]][solut->s[k_next]];

          cost_new = seq_get_C(solut,  0,  i_prev)                                                                  +   //       1st subseq
                  seq_get_W(solut,  j_next,  k)          * cost_concat_1 + seq_get_C(solut,  j_next,  k)             +   // concat 2nd subseq
                  seq_get_W(solut,  i,  j)               * cost_concat_2 + seq_get_C(solut,  i,  j)                  +   // concat 3rd subseq (reinserted seq)
                  seq_get_W(solut,  k_next,  data->dimen) * cost_concat_3 + seq_get_C(solut,  k_next,  data->dimen);       // concat 4th subseq
          
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
        update_subseq_info_matrix_b(solut, data, I < POS+1 ? I : POS+1);

        assert(abs(cost_best-solut->cost) < DBL_EPSILON);
        assert(feasible(solut->s, data->dimen));
        return true;
    }

    return false;
}


void RVND(tSolution * solut, tData * data) {

    int neighbd_list[] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    int nl_size = 5;
    int index;
    int neighbd;
    _Bool improve_flag;

    while (nl_size > 0) {
        index = data->rnd[data->rnd_index++];
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";
        //
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

        assert(feasible(solut->s, data->dimen));
      //if (feasible(solut->s, data->dimen+1)) {
      //    printf("qebrad\n");
      //    exit(0);
      //}
        //std::cout << (improve_flag ? "True" : "False") << std::endl;
        if (improve_flag) {
            neighbd_list[0] = SWAP;
            neighbd_list[1] = TWO_OPT;
            neighbd_list[2] = REINSERTION;
            neighbd_list[3] = OR_OPT_2;
            neighbd_list[4] = OR_OPT_3;
            nl_size = 5;

            //print_s(solut->s, info->dimen+1);
        } else {
            //print_s(neighbd_list, nl_size);
            memmove(neighbd_list + index, neighbd_list + index+1, sizeof(int)*(nl_size-index-1));
            nl_size--;
        }


    }

}

void perturb(tSolution * solut_crnt, tSolution * solut_partial, tData * data) {
    int s[data->dimen+1];
    memmove(s, solut_partial->s, sizeof(int)*(data->dimen+1));

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = floor((data->dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;
    //std::cout << "perturbing\n";
    //print_s(s);
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


    memcpy(solut_crnt->s, s, sizeof(int)*(data->dimen+1));
}


void GILS_RVND(int Imax, int Iils, tData * data) {

    tSolution solut_partial = Solution_init(*data);
    tSolution solut_crnt = Solution_init(*data);
    tSolution solut_best = Solution_init(*data);

#ifdef MATRIX
    printf("MATRIX\n");
#elif defined(FLAT)
    printf("FLAT\n");
#endif
    printf("struct tSolution's size: %zu bytes\n", solut_crnt.size);
    printf("struct tSolution's size: %.4f MB\n", solut_crnt.MBsize);


    for(int i = 0; i < Imax; ++i){
        int aux = data->rnd[data->rnd_index++];

        Real alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	


        construct(solut_crnt.s, alpha, data);
        assert(feasible(solut_crnt.s, data->dimen));
        print_s(solut_crnt.s, data->dimen+1);
        update_subseq_info_matrix(&solut_crnt, data);

        Solution_cpy(&solut_crnt, &solut_partial, data);
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);	

        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(&solut_crnt, data);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                Solution_cpy(&solut_crnt, &solut_partial, data);
                iterILS = 0;
            }

            perturb(&solut_crnt, &solut_partial, data);
            update_subseq_info_matrix(&solut_crnt, data);

            assert(feasible(solut_crnt.s, data->dimen));
            iterILS++;
        }

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(&solut_partial, &solut_best, data);
        }

        printf("\tCurrent best cost: %.2lf\n", solut_best.cost);

        printf("SOLUCAO: ");
        for(int i = 0; i < data->dimen+1; i++){
            printf("%d ", solut_best.s[i]);
        }
        printf("\n");

    }
    printf("COST: %.2lf\n", solut_best.cost);

    Solution_free(&solut_best);
    Solution_free(&solut_crnt);
    Solution_free(&solut_partial);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tData data = {0};

    int * rnd;

    data.dimen = loadData(&data.cost, &rnd);
    data.rnd = rnd;
    data.rnd_index = 0;

    Iils = data.dimen < 100 ? data.dimen : 100;
    time_t start = clock();
    GILS_RVND(Imax, Iils, &data);

    double res = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("TIME: %.6lf\n", res);

    tData_free(&data);

    return 0;
}

