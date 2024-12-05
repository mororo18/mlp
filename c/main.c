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

#define MATRIX

typedef unsigned uint;

typedef struct tData {
    double ** cost;
    int dimen;
    int * rnd;
    uint rnd_index;
} tData;

typedef struct tInfo {
    uint T;
    uint C;
    uint W;
} tInfo;

typedef struct tSolution {
    double *** seq;
    double cost;
    int * s;
    //int s_size;
} tSolution;

tSolution Solution_init(tInfo info, tData data) {
    tSolution solut;
    solut.s = (int*) calloc(data.dimen+1, sizeof(int));
    //solut.s_size = data.dimen+1;

    solut.seq = (double ***) calloc(data.dimen+1, sizeof(double **));
    for (int i = 0; i < data.dimen+1; i++) {
        solut.seq[i] = (double **) calloc(data.dimen+1, sizeof(double *));
        for (int j = 0; j < data.dimen+1; j++) {
            solut.seq[i][j] = (double *) calloc(3, sizeof(double));
        }
    }

    /*
    solut.seq = new tSeqInfo * [data.dimen+1];
    for (int i = 0; i < data.dimen+1; i++) {
        solut.seq[i] = new tSeqInfo [data.dimen+1];
        memset(solut.seq[i], 0.0, (data.dimen+1)*sizeof(tSeqInfo));
    }
    */

    solut.cost = DBL_MAX;

    return solut;
}

void Solution_cpy(tSolution * src, tSolution * tgt, const tData * data) {

    memcpy(tgt->s, src->s, sizeof(int)*(data->dimen+1));
    tgt->cost = src->cost;

    /*
    for (int i = 0; i < data.dimen+1; i++) {
        for (int j = 0; j < data.dimen+1; j++) {
            //memcpy(tgt.seq[i][j], src.seq[i][j], 3 * sizeof(double));
            std::copy(src.seq[i][j], src.seq[i][j] + 3, tgt.seq[i][j]);
        }
    }
    */

}

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


int partition(int *arr, int left, int right, tData * data, int r) {
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

void construct(int ** ret, const double alpha, tData * data){

    int s[data->dimen+1];

    memset(s, 0, (data->dimen+1)*sizeof(int));

    int cL[data->dimen-1];
    int cL_size = data->dimen-1;
    
    for (int i = 1; i < data->dimen; ++i) {
        cL[i-1] = i;
    }

    int r = 0;
    for (int j = 1; j < data->dimen; ++j) {
        sort(cL, cL_size, r, data);

        /**/
        int range = ceil(cL_size * alpha);
        int index = range > 0 ? rand() % range : 0;
        /**/

        //std::cout << data.rnd[data.rnd_index]<< std::endl;
        index = data->rnd[data->rnd_index++];
        int c = cL[index];
        s[j] = c;
        //print_s(cL);
        r = c;
        memmove(cL+index, cL+index+1, sizeof(int)*(cL_size-index));
        cL_size--;
    }

    memcpy(*ret, s, sizeof(int)*(data->dimen+1));
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

void update_subseq_info_matrix(tSolution * solut, const tInfo * info, const tData * data){
    int i, j, j_prev, k;
    int from = 0;
    char t;
    for (i = 0; i < data->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut->seq[i][i][info->T] = 0.0;
        solut->seq[i][i][info->C] = 0.0;
        solut->seq[i][i][info->W] = (double) !(i == 0);

        for (j = i+1; j < data->dimen+1; j++) {
            j_prev = j-1;
            
            solut->seq[i][j][info->T] = data->cost[solut->s[j_prev]][solut->s[j]] + solut->seq[i][j_prev][info->T];
            solut->seq[i][j][info->C] = solut->seq[i][j][info->T] + solut->seq[i][j_prev][info->C];
            solut->seq[i][j][info->W] = j + k;

        }
        from += t;
    }

    solut->cost = solut->seq[0][data->dimen][info->C];
}

void update_subseq_info_matrix_b(tSolution * solut, const tInfo * info, const tData * data, int index){
    int i, j, j_prev, k;
    int from = index;
    char t;
    for (i = 0; i < data->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut->seq[i][i][info->T] = 0.0;
        solut->seq[i][i][info->C] = 0.0;
        solut->seq[i][i][info->W] = (double) !(i == 0);

        for (j = i+1; j < data->dimen+1; j++) {
            j_prev = j-1;
            
            solut->seq[i][j][info->T] = data->cost[solut->s[j_prev]][solut->s[j]] + solut->seq[i][j_prev][info->T];
            solut->seq[i][j][info->C] = solut->seq[i][j][info->T] + solut->seq[i][j_prev][info->C];
            solut->seq[i][j][info->W] = j + k;

        }
        from += t;
    }

    solut->cost = solut->seq[0][data->dimen][info->C];
}

char search_swap(tSolution * solut, const tInfo * info, const tData * data) {
    double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    double cost_best = DBL_MAX;
    int i, j, j_prev, j_next, i_prev, i_next;
    int I;
    int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut->seq[0][i_prev][info->T] + data->cost[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + solut->seq[i][i_next][info->T]  + data->cost[solut->s[i]][solut->s[i_next+1]];

        cost_new = solut->seq[0][i_prev][info->C]                                                    +           //       1st subseq
        solut->seq[i][i_next][info->W]               * (cost_concat_1) + data->cost[solut->s[i_next]][solut->s[i]]  +           // concat 2nd subseq
        solut->seq[i_next+1][data->dimen][info->W]   * (cost_concat_2) + solut->seq[i_next+1][data->dimen][info->C];   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data->dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut->seq[0][i_prev][info->T]       + data->cost[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1                           + data->cost[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + solut->seq[i_next][j_prev][info->T]  + data->cost[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3                           + data->cost[solut->s[i]][solut->s[j_next]];

            cost_new = solut->seq[0][i_prev][info->C]                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut->seq[i_next][j_prev][info->W]      * cost_concat_2 + solut->seq[i_next][j_prev][info->C] +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut->seq[j_next][data->dimen][info->W] * cost_concat_4 + solut->seq[j_next][data->dimen][info->C];   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        swap(solut->s, I, J);
        update_subseq_info_matrix_b(solut, info, data, I);
        return TRUE;
    }

    return FALSE;
}

char search_two_opt(tSolution * solut, const tInfo * info, const tData * data) {
    double cost_new, 
        cost_concat_1, cost_concat_2;
    double cost_best = DBL_MAX;// cost_l1, cost_l2;
    double rev_seq_cost;
    int i, j, i_prev, j_next;
    int I;
    int J;

    for (i = 1; i < data->dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut->seq[i][i+1][info->T];
        for (j = i + 2; j < data->dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += data->cost[solut->s[j-1]][solut->s[j]] * (solut->seq[i][j][info->W]-1.0);

          cost_concat_1 =                 solut->seq[0][i_prev][info->T]   + data->cost[solut->s[j]][solut->s[i_prev]];
          cost_concat_2 = cost_concat_1 + solut->seq[i][j][info->T]        + data->cost[solut->s[j_next]][solut->s[i]];

          cost_new = solut->seq[0][i_prev][info->C]                                                        +   //  1st subseq
              solut->seq[i][j][info->W]                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut->seq[j_next][data->dimen][info->W] * cost_concat_2 + solut->seq[j_next][data->dimen][info->C];      // concat 3rd subseq

            
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reverse(solut->s, I, J);
        update_subseq_info_matrix(solut, info, data);
        return TRUE;
    }

    return FALSE;
}

char search_reinsertion(tSolution * solut, const tInfo * info, const tData * data, const int opt) {
    double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
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

          cost_concat_1 =                 solut->seq[0][k][info->T]            + data->cost[solut->s[k]][solut->s[i]];
          cost_concat_2 = cost_concat_1 + solut->seq[i][j][info->T]            + data->cost[solut->s[j]][solut->s[k_next]];
          cost_concat_3 = cost_concat_2 + solut->seq[k_next][i_prev][info->T]  + data->cost[solut->s[i_prev]][solut->s[j_next]];

          cost_new = solut->seq[0][k][info->C]                                                                   +   //       1st subseq
              solut->seq[i][j][info->W]               * cost_concat_1 + solut->seq[i][j][info->C]                  +   //  concat 2nd subseq (reinserted seq)
              solut->seq[k_next][i_prev][info->W]     * cost_concat_2 + solut->seq[k_next][i_prev][info->C]        +   //  concat 3rd subseq
              solut->seq[j_next][data->dimen][info->W] * cost_concat_3 + solut->seq[j_next][data->dimen][info->C];       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data->dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut->seq[0][i_prev][info->T]  + data->cost[solut->s[i_prev]][solut->s[j_next]];
          cost_concat_2 = cost_concat_1 + solut->seq[j_next][k][info->T]  + data->cost[solut->s[k]][solut->s[i]];
          cost_concat_3 = cost_concat_2 + solut->seq[i][j][info->T]       + data->cost[solut->s[j]][solut->s[k_next]];

          cost_new = solut->seq[0][i_prev][info->C]                                                                  +   //       1st subseq
                  solut->seq[j_next][k][info->W]          * cost_concat_1 + solut->seq[j_next][k][info->C]             +   // concat 2nd subseq
                  solut->seq[i][j][info->W]               * cost_concat_2 + solut->seq[i][j][info->C]                  +   // concat 3rd subseq (reinserted seq)
                  solut->seq[k_next][data->dimen][info->W] * cost_concat_3 + solut->seq[k_next][data->dimen][info->C];       // concat 4th subseq
          
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
        update_subseq_info_matrix_b(solut, info, data, I < POS+1 ? I : POS+1);
        return TRUE;
    }

    return FALSE;
}


void RVND(tSolution * solut, tInfo * info, tData * data) {

    int neighbd_list[] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    int nl_size = 5;
    uint index;
    int neighbd;
    char improve_flag;

    //cout << "RVND" << endl;
    while (nl_size > 0) {
        //k++;

        index = rand() % nl_size;
        index = data->rnd[data->rnd_index++];
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";

        improve_flag = FALSE;

        switch(neighbd){
            case REINSERTION:
                //before();
                improve_flag = search_reinsertion(solut, info, data, REINSERTION);
                //after(REINSERTION);
                break;				
            case OR_OPT_2:
                //before();
                improve_flag = search_reinsertion(solut, info, data, OR_OPT_2);
                //after(OR_OPT2);
                break;				
            case OR_OPT_3:
                //before();
                improve_flag = search_reinsertion(solut, info, data, OR_OPT_3);
                //after(OR_OPT3);
                break;				
            case SWAP:
                //before();
                improve_flag = search_swap(solut, info, data);
                //after(SWAP);
                break;
            case TWO_OPT:
                //before();
                improve_flag = search_two_opt(solut, info, data);
                //after(TWO_OPT);
                break;				
        }
        //std::cout << (improve_flag ? "True" : "False") << std::endl;
        if (improve_flag) {
            neighbd_list[0] = SWAP;
            neighbd_list[1] = TWO_OPT;
            neighbd_list[2] = REINSERTION;
            neighbd_list[3] = OR_OPT_2;
            neighbd_list[4] = OR_OPT_3;
            nl_size = 5;
        } else {
            //std::cout << index << "  " << neighbd_list.size() << std::endl;
            //std::cout << solut.cost << std::endl;
            
            //std::cout << data.rnd_index << std::endl;
            memmove(neighbd_list + index, neighbd_list + index+1, sizeof(int)*(nl_size-index-1));
            nl_size--;
        }

        //std::cout << "cost  " << solut.cost << std::endl ;


    }

    //exit(0);
    //std::cout << k << " RVND iteracoes" << std::endl;
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
        /**/
        int max = (data->dimen+1) -2 -size_max;
        A_start = rand() % max + 1;
        A_end = A_start + rand() % (size_max - size_min + 1) + size_min;

        B_start = rand() % max + 1;
        B_end = B_start + rand() % (size_max - size_min + 1) + size_min;
        /**/



        //std::cout << "paa\n";

        //cout << data.rnd[data.rnd_index] << endl;
        A_start = data->rnd[data->rnd_index++];
        //cout << data.rnd[data.rnd_index] << endl;
        A_end = A_start + data->rnd[data->rnd_index++];
        //std::cout << "A start  " << A_start << std::endl;
        //std::cout << "A end  " << A_end << std::endl;

        //cout << data.rnd[data.rnd_index] << endl;
        B_start = data->rnd[data->rnd_index++];
        //cout << data.rnd[data.rnd_index] << endl;
        B_end = B_start + data->rnd[data->rnd_index++];
        //std::cout << "B start  " << B_start << std::endl;
        //std::cout << "B end  " << B_end << std::endl;
    }
    
    //cout << "A_end  " << A_end << endl << "B_end  " << B_end << endl;

    if (A_start < B_start) {
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    } else {
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    //print_s(s);
    //update_subseq_info_matrix(solut, info);

    memcpy(solut_crnt->s, s, sizeof(int)*(data->dimen+1));
}


void GILS_RVND(int Imax, int Iils, tInfo * info, tData * data) {

    tSolution solut_partial = Solution_init(*info, *data);
    tSolution solut_crnt = Solution_init(*info, *data);
    tSolution solut_best = Solution_init(*info, *data);

    for(int i = 0; i < Imax; ++i){
        /**/ int aux = (unsigned)rand() % TABLE_SZ;
        aux = data->rnd[data->rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	


        construct(&solut_crnt.s, alpha, data);
        //print_s(solut_crnt.s);
        update_subseq_info_matrix(&solut_crnt, info, data);

        //solut_partial = solut_crnt;
        Solution_cpy(&solut_crnt, &solut_partial, data);
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            RVND(&solut_crnt, info, data);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                Solution_cpy(&solut_crnt, &solut_partial, data);
                //solut_partial = solut_crnt;
                iterILS = 0;
            }

            perturb(&solut_crnt, &solut_partial, data);
            update_subseq_info_matrix(&solut_crnt, info, data);
            //exit(0);
            //std::cout << "ITER  " << iterILS << std::endl;
            iterILS++;
        }

        //update_subseq_info_matrix(solut_partial, info);

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(&solut_partial, &solut_best, data);
            //solut_best = solut_partial;
        }

        //after(7);

        //std::cout << "\tCurrent search cost: "<< cost_sl << std::endl;
        printf("\tCurrent best cost: %.2lf\n", solut_best.cost);
        //std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        //std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;
        //std::cout << k << "  Iteracoes " << std::endl;

        printf("SOLUCAO: ");
        for(int i = 0; i < data->dimen+1; i++){
            printf("%d ", solut_best.s[i]);
        }
        printf("\n");

    }
    //std::cout << "Dimension: " << dimension << std::endl;
    printf("COST: %.2lf\n", solut_best.cost);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tInfo info = {0};
    info.T = 0;
    info.W = 1;
    info.C = 2;

    int * rnd;

    tData data = {0};
    data.dimen = loadData(&data.cost, &rnd);
    data.rnd = rnd;
    //print_s(rnd);
    //printf("%d\n", rnd[10]);
    data.rnd_index = 0;

    tSolution solut = Solution_init(info, data);

    //exit(0);
    /*
    for (int i = 0; i <=n; i++) {
        for (int j = i+1; j <=n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */


    srand(clock());

    Iils = data.dimen < 100 ? data.dimen : 100;
    time_t start = clock();
    GILS_RVND(Imax, Iils, &info, &data);

    double res = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("TIME: %.6lf\n", res);

    /*
    flag = true;

    srand(duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count());
    readData(argc, argv, &dimension, &c);

    int ar[] = {100, dimension};
    
    Iils = ar[dimension < 100];

    auto t1 = high_resolution_clock::now();

    GILS_RVND(Imax, Iils);

    auto t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    double res = (double)duration / 10e2;
    std::cout << "TIME: " << res << std::endl;

    if(flag){
        std::cout << "Construction time: " << construct_t/10e5 << std::endl;
        std::cout << "Swap time: " << swap_t/10e5 << std::endl;
        std::cout << "two_opt time: " << two_opt_t/10e5 << std::endl;
        std::cout << "reinsertion time: " << reinsertion_t/10e5 << std::endl;
        std::cout << "or_opt2 time: " << opt2_t/10e5 << std::endl;
        std::cout << "or_opt3 time: " << opt3_t /10e5<< std::endl;
    }
    */
    return 0;
}

