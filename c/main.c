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
    uint rnd_index;
} tRnd;

typedef struct tSubseq {
    double C;
    double T;
    double W;
} tSubseq;


//  tSolution Solution_init(tInfo info) {
//      tSolution solut;
//      solut.s = (int*) calloc(info.dimen+1, sizeof(int));
//      //solut.s_size = info.dimen+1;

//      solut.seq = (double ***) calloc(info.dimen+1, sizeof(double **));
//      for (int i = 0; i < info.dimen+1; i++) {
//          solut.seq[i] = (double **) calloc(info.dimen+1, sizeof(double *));
//          for (int j = 0; j < info.dimen+1; j++) {
//              solut.seq[i][j] = (double *) calloc(3, sizeof(double));
//          }
//      }

//      /*
//      solut.seq = new tSeqInfo * [info.dimen+1];
//      for (int i = 0; i < info.dimen+1; i++) {
//          solut.seq[i] = new tSeqInfo [info.dimen+1];
//          memset(solut.seq[i], 0.0, (info.dimen+1)*sizeof(tSeqInfo));
//      }
//      */

//      solut.cost = DBL_MAX;

//      return solut;
//  }

//  void Solution_cpy(tSolution * src, tSolution * tgt, const tInfo * info) {

//      memcpy(tgt->s, src->s, sizeof(int)*(info->dimen+1));
//      tgt->cost = src->cost;

//      /*
//      for (int i = 0; i < info.dimen+1; i++) {
//          for (int j = 0; j < info.dimen+1; j++) {
//              //memcpy(tgt.seq[i][j], src.seq[i][j], 3 * sizeof(double));
//              std::copy(src.seq[i][j], src.seq[i][j] + 3, tgt.seq[i][j]);
//          }
//      }
//      */

//  }

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

void sort(int * arr, int arr_size, int r) {

    for (int i = 0; i < arr_size; i++) {
        for (int j = 0; j < arr_size-i-1; j++) {
            if (cost[r][arr[j]] > cost[r][arr[j+1]]) {
                int tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;
            }
        }
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

        /**/
        int range = ceil(cL_size * alpha);
        int index = range > 0 ? rand() % range : 0;
        /**/

        //std::cout << info.rnd[info.rnd_index]<< std::endl;
        index = rnd->rnd[rnd->rnd_index++];
        int c = cL[index];
        s[j] = c;
        //print_s(cL);
        r = c;
        memmove(cL+index, cL+index+1, sizeof(int)*(cL_size-index));
        cL_size--;
    }

    memcpy(ret, s, sizeof(int)*(dimen+1));
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

double subseq_load(int * s, tSubseq seq[][dimen+1]){
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

double subseq_load_b(int * s, tSubseq seq[][dimen+1], int index){
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

char search_swap(tSolution * solut, const tInfo * info) {
    double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    double cost_best = DBL_MAX;
    int i, j, j_prev, j_next, i_prev, i_next;
    int I;
    int J;

    for (i = 1; i < info->dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut->seq[0][i_prev][info->T] + info->cost[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + solut->seq[i][i_next][info->T]  + info->cost[solut->s[i]][solut->s[i_next+1]];

        cost_new = solut->seq[0][i_prev][info->C]                                                    +           //       1st subseq
        solut->seq[i][i_next][info->W]               * (cost_concat_1) + info->cost[solut->s[i_next]][solut->s[i]]  +           // concat 2nd subseq
        solut->seq[i_next+1][info->dimen][info->W]   * (cost_concat_2) + solut->seq[i_next+1][info->dimen][info->C];   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < info->dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut->seq[0][i_prev][info->T]       + info->cost[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1                           + info->cost[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + solut->seq[i_next][j_prev][info->T]  + info->cost[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3                           + info->cost[solut->s[i]][solut->s[j_next]];

            cost_new = solut->seq[0][i_prev][info->C]                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut->seq[i_next][j_prev][info->W]      * cost_concat_2 + solut->seq[i_next][j_prev][info->C] +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut->seq[j_next][info->dimen][info->W] * cost_concat_4 + solut->seq[j_next][info->dimen][info->C];   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        swap(solut->s, I, J);
        subseq_load_b(solut, info, I);
        return TRUE;
    }

    return FALSE;
}

char search_two_opt(tSolution * solut, const tInfo * info) {
    double cost_new, 
        cost_concat_1, cost_concat_2;
    double cost_best = DBL_MAX;// cost_l1, cost_l2;
    double rev_seq_cost;
    int i, j, i_prev, j_next;
    int I;
    int J;

    for (i = 1; i < info->dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut->seq[i][i+1][info->T];
        for (j = i + 2; j < info->dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += info->cost[solut->s[j-1]][solut->s[j]] * (solut->seq[i][j][info->W]-1.0);

          cost_concat_1 =                 solut->seq[0][i_prev][info->T]   + info->cost[solut->s[j]][solut->s[i_prev]];
          cost_concat_2 = cost_concat_1 + solut->seq[i][j][info->T]        + info->cost[solut->s[j_next]][solut->s[i]];

          cost_new = solut->seq[0][i_prev][info->C]                                                        +   //  1st subseq
              solut->seq[i][j][info->W]                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut->seq[j_next][info->dimen][info->W] * cost_concat_2 + solut->seq[j_next][info->dimen][info->C];      // concat 3rd subseq

            
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost - DBL_EPSILON) {
        reverse(solut->s, I, J);
        subseq_load(solut, info);
        return TRUE;
    }

    return FALSE;
}

char search_reinsertion(tSolution * solut, const tInfo * info, const int opt) {
    double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    int i, j, k, k_next, i_prev, j_next;
    int I;
    int J;
    int POS;

    for (i = 1, j = opt +i-1; i < info->dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 solut->seq[0][k][info->T]            + info->cost[solut->s[k]][solut->s[i]];
          cost_concat_2 = cost_concat_1 + solut->seq[i][j][info->T]            + info->cost[solut->s[j]][solut->s[k_next]];
          cost_concat_3 = cost_concat_2 + solut->seq[k_next][i_prev][info->T]  + info->cost[solut->s[i_prev]][solut->s[j_next]];

          cost_new = solut->seq[0][k][info->C]                                                                   +   //       1st subseq
              solut->seq[i][j][info->W]               * cost_concat_1 + solut->seq[i][j][info->C]                  +   //  concat 2nd subseq (reinserted seq)
              solut->seq[k_next][i_prev][info->W]     * cost_concat_2 + solut->seq[k_next][i_prev][info->C]        +   //  concat 3rd subseq
              solut->seq[j_next][info->dimen][info->W] * cost_concat_3 + solut->seq[j_next][info->dimen][info->C];       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < info->dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut->seq[0][i_prev][info->T]  + info->cost[solut->s[i_prev]][solut->s[j_next]];
          cost_concat_2 = cost_concat_1 + solut->seq[j_next][k][info->T]  + info->cost[solut->s[k]][solut->s[i]];
          cost_concat_3 = cost_concat_2 + solut->seq[i][j][info->T]       + info->cost[solut->s[j]][solut->s[k_next]];

          cost_new = solut->seq[0][i_prev][info->C]                                                                  +   //       1st subseq
                  solut->seq[j_next][k][info->W]          * cost_concat_1 + solut->seq[j_next][k][info->C]             +   // concat 2nd subseq
                  solut->seq[i][j][info->W]               * cost_concat_2 + solut->seq[i][j][info->C]                  +   // concat 3rd subseq (reinserted seq)
                  solut->seq[k_next][info->dimen][info->W] * cost_concat_3 + solut->seq[k_next][info->dimen][info->C];       // concat 4th subseq
          
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
        subseq_load_b(solut, info, I < POS+1 ? I : POS+1);
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

    //cout << "RVND" << endl;
    while (nl_size > 0) {
        //k++;

        index = rand() % nl_size;
        index = rnd->rnd[rnd->rnd_index++];
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";

        improve_flag = FALSE;

        switch(neighbd){
            case REINSERTION:
                //before();
                improve_flag = search_reinsertion(solut, seq, REINSERTION);
                //after(REINSERTION);
                break;				
            case OR_OPT_2:
                //before();
                improve_flag = search_reinsertion(solut, seq, OR_OPT_2);
                //after(OR_OPT2);
                break;				
            case OR_OPT_3:
                //before();
                improve_flag = search_reinsertion(solut, seq, OR_OPT_3);
                //after(OR_OPT3);
                break;				
            case SWAP:
                //before();
                improve_flag = search_swap(solut, seq);
                //after(SWAP);
                break;
            case TWO_OPT:
                //before();
                improve_flag = search_two_opt(solut, seq);
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
            
            //std::cout << info.rnd_index << std::endl;
            memmove(neighbd_list + index, neighbd_list + index+1, sizeof(int)*(nl_size-index-1));
            nl_size--;
        }

        //std::cout << "cost  " << solut.cost << std::endl ;


    }

    //exit(0);
    //std::cout << k << " RVND iteracoes" << std::endl;
}

void perturb(int * s_crnt, int * s_partial, tRnd * rnd) {
    int s[dimen+1];
    memmove(s, s_partial, sizeof(int)*(dimen+1));

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = floor((dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;
    //std::cout << "perturbing\n";
    //print_s(s);
    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        /**/
        int max = (dimen+1) -2 -size_max;
        A_start = rand() % max + 1;
        A_end = A_start + rand() % (size_max - size_min + 1) + size_min;

        B_start = rand() % max + 1;
        B_end = B_start + rand() % (size_max - size_min + 1) + size_min;
        /**/



        //std::cout << "paa\n";

        //cout << info.rnd[info.rnd_index] << endl;
        A_start = rnd->rnd[rnd->rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        A_end = A_start + rnd->rnd[rnd->rnd_index++];
        //std::cout << "A start  " << A_start << std::endl;
        //std::cout << "A end  " << A_end << std::endl;

        //cout << info.rnd[info.rnd_index] << endl;
        B_start = rnd->rnd[rnd->rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        B_end = B_start + rnd->rnd[rnd->rnd_index++];
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
    //subseq_load(solut, info);

    memcpy(s_crnt, s, sizeof(int)*(dimen+1));
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
        /**/ int aux = (unsigned)rand() % TABLE_SZ;
        aux = rnd->rnd[rnd->rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	

        construct(s_crnt, alpha, rnd);
        //print_s(solut_crnt.s);
        cost_crnt = subseq_load(s_crnt, seq);
        cost_partial = cost_crnt;

        //solut_partial = solut_crnt;
        memcpy(s_partial, s_crnt, sizeof(int)*(dimen+1));
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", cost_partial);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            RVND(s_crnt, seq);
            if(cost_crnt < cost_partial - DBL_EPSILON){
                Solution_cpy(&solut_crnt, &solut_partial, info);
                //solut_partial = solut_crnt;
                iterILS = 0;
            }

            perturb(&solut_crnt, &solut_partial, info);
            subseq_load(s_crnt, seq);
            //exit(0);
            //std::cout << "ITER  " << iterILS << std::endl;
            iterILS++;
        }

        //subseq_load(solut_partial, info);

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(&solut_partial, &solut_best, info);
            //solut_best = solut_partial;
        }

        //after(7);

        //std::cout << "\tCurrent search cost: "<< cost_sl << std::endl;
        printf("\tCurrent best cost: %.2lf\n", solut_best.cost);
        //std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        //std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;
        //std::cout << k << "  Iteracoes " << std::endl;

        printf("SOLUCAO: ");
        for(int i = 0; i < info->dimen+1; i++){
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

    tRnd rnd = {0}

    int * rnd;

    dimen = loadData(&cost, &rnd);
    rnd.rnd = rnd;
    rnd.rnd_index = 0;
    //print_s(rnd);
    //printf("%d\n", rnd[10]);

    srand(clock());

    Iils = dimen < 100 ? dimen : 100;
    time_t start = clock();
    GILS_RVND(Imax, Iils, &rnd);

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

