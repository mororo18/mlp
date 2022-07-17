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

typedef struct tInfo {
    double ** cost;
    int dimen;
    uint T;
    uint C;
    uint W;
    int * rnd;
    uint rnd_index;
} tInfo;

typedef struct tSeqInfo {
    double T, C, W;
} tSeqInfo;

typedef tSeqInfo tSeq;
typedef tSeq * tSeq_;
typedef tSeq_ * tSeq__;

typedef struct tSolution {
    tSeq__ seq;
    //double *** seq;
    double cost;
    int * s;
    //int s_size;
} tSolution;

inline void set_C(tSeq__ seq, int i, int j, double value) {
    seq[i][j].C = value;
}

inline void set_T(tSeq__ seq, int i, int j, double value) {
    seq[i][j].T = value;
}

inline void set_W(tSeq__ seq, int i, int j, double value) {
    seq[i][j].W = value;
}

inline double get_C(tSeq__ seq, int i, int j) {
    return seq[i][j].C;
}

inline double get_T(tSeq__ seq, int i, int j) {
    return seq[i][j].T;
}

inline double get_W(tSeq__ seq, int i, int j) {
    return seq[i][j].W;
}

tSolution Solution_init(tInfo info) {
    tSolution solut;
    solut.s = (int*) calloc(info.dimen+1, sizeof(int));
    //solut.s_size = info.dimen+1;

  //solut.seq = (double ***) calloc(info.dimen+1, sizeof(double **));
  //for (int i = 0; i < info.dimen+1; i++) {
  //    solut.seq[i] = (double **) calloc(info.dimen+1, sizeof(double *));
  //    for (int j = 0; j < info.dimen+1; j++) {
  //        solut.seq[i][j] = (double *) calloc(3, sizeof(double));
  //    }
  //}

    solut.seq = (tSeq__) calloc(info.dimen+1, sizeof(tSeq_));
    for (int i = 0; i < info.dimen+1; i++) {
        solut.seq[i] = (tSeq_) calloc(info.dimen+1, sizeof(tSeq));
      //for (int j = 0; j < info.dimen+1; j++) {
      //    solut.seq[i][j] = (double *) calloc(3, sizeof(double));
      //}
    }

    solut.cost = DBL_MAX;

    return solut;
}

void Solution_cpy(tSolution * src, tSolution * tgt, const tInfo * info) {

    memcpy(tgt->s, src->s, sizeof(int)*(info->dimen+1));
    tgt->cost = src->cost;

    /*
    for (int i = 0; i < info.dimen+1; i++) {
        for (int j = 0; j < info.dimen+1; j++) {
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
        printf("%d ", s[i]);
    printf("\n");
}

void sort(int * arr, int arr_size, int r, tInfo * info) {

    for (int i = 0; i < arr_size; i++) {
        for (int j = 0; j < arr_size-i-1; j++) {
            if (info->cost[r][arr[j]] > info->cost[r][arr[j+1]]) {
                int tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;
            }
        }
    }
}

void construct(int * ret, const double alpha, tInfo * info){

    int s[info->dimen+1];

    memset(s, 0, (info->dimen+1)*sizeof(int));

    int cL[info->dimen-1];
    int cL_size = info->dimen-1;
    
    for (int i = 1; i < info->dimen; ++i) {
        cL[i-1] = i;
    }

    int r = 0;
    for (int j = 1; j < info->dimen; ++j) {
        sort(cL, cL_size, r, info);

        /**/
        int range = ceil(cL_size * alpha);
        int index = range > 0 ? rand() % range : 0;
        /**/

        //std::cout << info.rnd[info.rnd_index]<< std::endl;
        index = info->rnd[info->rnd_index++];
        int c = cL[index];
        s[j] = c;
        //print_s(cL);
        r = c;
        memmove(cL+index, cL+index+1, sizeof(int)*(cL_size-index));
        cL_size--;
    }

    memcpy(ret, s, sizeof(int)*(info->dimen+1));
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

void subseq_load(tSolution * solut, const tInfo * info){
    int i, j, j_prev, k;
    int from = 0;
    char t;
    for (i = 0; i < info->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        set_T(solut->seq, i, i, 0.0);
        set_C(solut->seq, i, i, 0.0);
        set_W(solut->seq, i, i, (double)!(i == 0));

        for (j = i+1; j < info->dimen+1; j++) {
            j_prev = j-1;
            
            double T = info->cost[solut->s[j_prev]][solut->s[j]] + get_T(solut->seq, i, j_prev);
            set_T(solut->seq, i, j, T);

            double C = get_T(solut->seq, i, j) + get_C(solut->seq, i, j_prev);
            set_C(solut->seq, i, j, C);

            double W = j + k;
            set_W(solut->seq, i, j, W);

        }
        from += t;
    }

    solut->cost = get_C(solut->seq, 0, info->dimen);
}

void subseq_load_b(tSolution * solut, const tInfo * info, int index){
    int i, j, j_prev, k;
    int from = index;
    char t;
    for (i = 0; i < info->dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        set_T(solut->seq, i, i, 0.0);
        set_C(solut->seq, i, i, 0.0);
        set_W(solut->seq, i, i, (double)!(i == 0));

        for (j = i+1; j < info->dimen+1; j++) {
            j_prev = j-1;
            
            double T = info->cost[solut->s[j_prev]][solut->s[j]] + get_T(solut->seq, i, j_prev);
            set_T(solut->seq, i, j, T);

            double C = get_T(solut->seq, i, j) + get_C(solut->seq, i, j_prev);
            set_C(solut->seq, i, j, C);

            double W = j + k;
            set_W(solut->seq, i, j, W);


        }
        from += t;
    }

    solut->cost = get_C(solut->seq, 0, info->dimen);
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
        cost_concat_1 =                 get_T(solut->seq, 0, i_prev) + info->cost[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + get_T(solut->seq, i, i_next)  + info->cost[solut->s[i]][solut->s[i_next+1]];

        cost_new = get_C(solut->seq, 0, i_prev)                                                    +           //       1st subseq
        get_W(solut->seq, i, i_next)               * (cost_concat_1) + info->cost[solut->s[i_next]][solut->s[i]]  +           // concat 2nd subseq
        get_W(solut->seq, i_next+1, info->dimen)   * (cost_concat_2) + get_C(solut->seq, i_next+1, info->dimen);   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < info->dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 get_T(solut->seq, 0, i_prev)       + info->cost[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1                           + info->cost[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + get_T(solut->seq, i_next, j_prev)  + info->cost[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3                           + info->cost[solut->s[i]][solut->s[j_next]];

            cost_new = get_C(solut->seq, 0, i_prev)                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            get_W(solut->seq, i_next, j_prev)      * cost_concat_2 + get_C(solut->seq, i_next, j_prev) +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            get_W(solut->seq, j_next, info->dimen) * cost_concat_4 + get_C(solut->seq, j_next, info->dimen);   // concat 5th subseq

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

        rev_seq_cost = get_T(solut->seq, i, i+1);
        for (j = i + 2; j < info->dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += info->cost[solut->s[j-1]][solut->s[j]] * (get_W(solut->seq, i, j)-1.0);

          cost_concat_1 =                 get_T(solut->seq, 0, i_prev)   + info->cost[solut->s[j]][solut->s[i_prev]];
          cost_concat_2 = cost_concat_1 + get_T(solut->seq, i, j)        + info->cost[solut->s[j_next]][solut->s[i]];

          cost_new = get_C(solut->seq, 0, i_prev)                                                        +   //  1st subseq
              get_W(solut->seq, i, j)                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              get_W(solut->seq, j_next, info->dimen) * cost_concat_2 + get_C(solut->seq, j_next, info->dimen);      // concat 3rd subseq

            
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

          cost_concat_1 =                 get_T(solut->seq, 0, k)            + info->cost[solut->s[k]][solut->s[i]];
          cost_concat_2 = cost_concat_1 + get_T(solut->seq, i, j)            + info->cost[solut->s[j]][solut->s[k_next]];
          cost_concat_3 = cost_concat_2 + get_T(solut->seq, k_next, i_prev)  + info->cost[solut->s[i_prev]][solut->s[j_next]];

          cost_new = get_C(solut->seq, 0, k)                                                                   +   //       1st subseq
              get_W(solut->seq, i, j)               * cost_concat_1 + get_C(solut->seq, i, j)                  +   //  concat 2nd subseq (reinserted seq)
              get_W(solut->seq, k_next, i_prev)     * cost_concat_2 + get_C(solut->seq, k_next, i_prev)        +   //  concat 3rd subseq
              get_W(solut->seq, j_next, info->dimen) * cost_concat_3 + get_C(solut->seq, j_next, info->dimen);       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < info->dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 get_T(solut->seq, 0, i_prev)  + info->cost[solut->s[i_prev]][solut->s[j_next]];
          cost_concat_2 = cost_concat_1 + get_T(solut->seq, j_next, k)  + info->cost[solut->s[k]][solut->s[i]];
          cost_concat_3 = cost_concat_2 + get_T(solut->seq, i, j)       + info->cost[solut->s[j]][solut->s[k_next]];

          cost_new = get_C(solut->seq, 0, i_prev)                                                                  +   //       1st subseq
                  get_W(solut->seq, j_next, k)          * cost_concat_1 + get_C(solut->seq, j_next, k)             +   // concat 2nd subseq
                  get_W(solut->seq, i, j)               * cost_concat_2 + get_C(solut->seq, i, j)                  +   // concat 3rd subseq (reinserted seq)
                  get_W(solut->seq, k_next, info->dimen) * cost_concat_3 + get_C(solut->seq, k_next, info->dimen);       // concat 4th subseq
          
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
        if (cost_best == 21575.0) {
            puts("pa");
        }
        return TRUE;
    }

    return FALSE;
}


char feasible(int * s, int sz) {
    int count[sz]; 

    memset(count, 0, sizeof(int)*sz);

    for (int i = 0; i < sz; i++)
        count[s[i]] += 1;

    for (int i = 0; i < sz; i++)
        if (count[i] != 1) 
            return FALSE;

    return TRUE;
}

void RVND(tSolution * solut, tInfo * info) {

    int neighbd_list[] = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    int nl_size = 5;
    uint index;
    int neighbd;
    char improve_flag;

    //printf("RVND\n");
    while (nl_size > 0) {
        //k++;
        //printf("\t%.0lf\n", solut->cost);
        index = rand() % nl_size;
        index = info->rnd[info->rnd_index++];
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";

        improve_flag = FALSE;

        switch(neighbd){
            case REINSERTION:
                //before();
                improve_flag = search_reinsertion(solut, info, REINSERTION);
                //printf("REINSERTION");
                break;				
            case OR_OPT_2:
                //before();
                improve_flag = search_reinsertion(solut, info, OR_OPT_2);
                //after();
                //printf("OR_OPT2\t");
                break;				
            case OR_OPT_3:
                //before();
                improve_flag = search_reinsertion(solut, info, OR_OPT_3);
                //after(OR_OPT3);
                //printf("OR_OPT3\t");
                break;				
            case SWAP:
                //before();
                improve_flag = search_swap(solut, info);
                //after(SWAP);
                //printf("SWAP\t");
                break;
            case TWO_OPT:
                //before();
                improve_flag = search_two_opt(solut, info);
                //after(TWO_OPT);
                //printf("TWO_OPT\t");
                break;				
        }

        if (feasible(solut->s, info->dimen+1)) {
            printf("qebrad\n");
            exit(0);
        }
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
            //std::cout << index << "  " << neighbd_list.size() << std::endl;
            //std::cout << solut.cost << std::endl;
            
            //std::cout << info.rnd_index << std::endl;
            //print_s(neighbd_list, nl_size);
            memmove(neighbd_list + index, neighbd_list + index+1, sizeof(int)*(nl_size-index-1));
            nl_size--;
        }

        //std::cout << "cost  " << solut.cost << std::endl ;


    }

    //exit(0);
    //std::cout << k << " RVND iteracoes" << std::endl;
}

void perturb(tSolution * solut_crnt, tSolution * solut_partial, tInfo * info) {
    int s[info->dimen+1];
    memmove(s, solut_partial->s, sizeof(int)*(info->dimen+1));

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = floor((info->dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;
    //std::cout << "perturbing\n";
    //print_s(s);
    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        /**/
        int max = (info->dimen+1) -2 -size_max;
        A_start = rand() % max + 1;
        A_end = A_start + rand() % (size_max - size_min + 1) + size_min;

        B_start = rand() % max + 1;
        B_end = B_start + rand() % (size_max - size_min + 1) + size_min;
        /**/



        //std::cout << "paa\n";

        //cout << info.rnd[info.rnd_index] << endl;
        A_start = info->rnd[info->rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        A_end = A_start + info->rnd[info->rnd_index++];
        //std::cout << "A start  " << A_start << std::endl;
        //std::cout << "A end  " << A_end << std::endl;

        //cout << info.rnd[info.rnd_index] << endl;
        B_start = info->rnd[info->rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        B_end = B_start + info->rnd[info->rnd_index++];
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

    memcpy(solut_crnt->s, s, sizeof(int)*(info->dimen+1));
}


void GILS_RVND(int Imax, int Iils, tInfo * info) {

    tSolution solut_partial = Solution_init(*info);
    tSolution solut_crnt = Solution_init(*info);
    tSolution solut_best = Solution_init(*info);

    for(int i = 0; i < Imax; ++i){
        /**/ int aux = (unsigned)rand() % TABLE_SZ;
        aux = info->rnd[info->rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	


        construct(solut_crnt.s, alpha, info);
        print_s(solut_crnt.s, info->dimen+1);
        subseq_load(&solut_crnt, info);

        //solut_partial = solut_crnt;
        Solution_cpy(&solut_crnt, &solut_partial, info);
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            RVND(&solut_crnt, info);
            //printf("%.2lf\n", solut_crnt.cost);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                Solution_cpy(&solut_crnt, &solut_partial, info);
                //solut_partial = solut_crnt;
                iterILS = 0;
            }

            perturb(&solut_crnt, &solut_partial, info);
            subseq_load(&solut_crnt, info);
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

    tInfo info = {0};
    info.T = 0;
    info.W = 1;
    info.C = 2;

    int * rnd;

    info.dimen = loadData(&info.cost, &rnd);
    info.rnd = rnd;
    //print_s(rnd);
    //printf("%d\n", rnd[10]);
    info.rnd_index = 0;

    tSolution solut = Solution_init(info);

    //exit(0);
    for (int i = 0; i < info.dimen; i++) {
        for (int j = 0; j < info.dimen; j++) {
            printf("%.0lf ",info.cost[i][j]);
        }
        puts("");
    }


    srand(clock());

    Iils = info.dimen < 100 ? info.dimen : 100;
    time_t start = clock();
    GILS_RVND(Imax, Iils, &info);

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

