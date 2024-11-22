#include <iostream>
#include <cstdint>
#include <cstring>
#include <cfloat>
#include <new>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "readData.h"
#include "Data.hpp"

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

#define FLAT

using std::chrono::high_resolution_clock;

typedef unsigned uint;

// 3D array

int dimen;
double ** cost;

typedef struct tRnd {
    vector<int> rnd;
    uint rnd_index;
} tRnd;

typedef struct seqStruct {
    double T;
    double C;
    double W;
} seqStruct;

typedef std::vector<std::vector<seqStruct>> tSubseq;

double R_table(int i){
    static const double table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

void print_s(std::vector<int> s) {

    for (int i = 0; i < s.size(); i++)
        std::cout << s[i]+1 << " ";
    std::cout << std::endl;
}

int partition(std::vector<int>& arr, int left, int right, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (cost[r][arr[j]] < cost[r][pivot]) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

void quicksort(std::vector<int>& arr, int left, int right, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, r);
        quicksort(arr, left, pivot - 1, r);
        quicksort(arr, pivot + 1, right, r);
    }
}

void sort(std::vector<int>& arr, int r) {
    quicksort(arr, 0, arr.size() - 1, r);
}

std::vector<int> construct(const double alpha, tRnd & rnd){

    std::vector<int> s = {0};
    s.reserve(dimen+1);
    std::vector<int> cL(dimen-1);
    for(int i = 1; i < dimen; ++i){
        cL[i-1] = i;
    }

    int r = 0;
    while (!cL.empty()) {
        sort(cL, r);

        /**/
        int range = std::ceil(cL.size() * alpha);
        int index = range > 0 ? rand() % range : 0;
        /**/

        //std::cout << info.rnd[info.rnd_index]<< std::endl;
        index = rnd.rnd[rnd.rnd_index++];
        int c = cL[index];
        s.push_back(c);
        //print_s(cL);
        r = c;
        cL.erase(cL.begin() + index);
    }

    s.push_back(0);

    return s;
}	

inline void swap(std::vector<int> &vec, int i, int j){
    std::iter_swap(vec.begin() + i, vec.begin() + j);
}

inline void reverse(std::vector<int> &vec, int i, int j){
    std::reverse(vec.begin() + i, vec.begin() + j+1);
}

inline void reinsert(std::vector<int> &vec, int i, int j, int pos){
    std::vector<int> seq (vec.begin() + i, vec.begin() +j+1);
    if(pos < i){
        vec.erase(vec.begin() + i, vec.begin() + j+1);
        vec.insert(vec.begin() + pos, seq.begin(), seq.end());
    }else{
        vec.insert(vec.begin() + pos, seq.begin(), seq.end());
        vec.erase(vec.begin() + i, vec.begin() + j+1);
    }

}

inline void subseq_load(vector<int> & s, tSubseq & seq, int index = 0){
    alignas(INT_SZ) int i, j, j_prev, k;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for (i = 0; i < dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        seq[i][ i].T = 0.0;
        seq[i][ i].C = 0.0;
        seq[i][ i].W = (double) !(i == 0);

        for (j = i+1; j < dimen+1; j++) {
            j_prev = j-1;
            
            seq[i][ j].T = cost[s[j_prev]][s[j]] + seq[i][j_prev].T;
            seq[i][j].C = seq[i][j].T + seq[i][j_prev].C;
            seq[i][ j].W = j + k;
        }
        from += t;
    }
}

bool improve = false;

void search_swap(vector<int> & s, tSubseq & seq) {
//inline bool search_swap(vector<int> & s, tSubseq & seq) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    alignas(DBL_SZ) double cost_best = DBL_MAX;
    alignas(INT_SZ) int i, j, j_prev, j_next, i_prev, i_next;
    //alignas(INT_SZ) int dim = dimension - 2;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 seq[0][ i_prev].T + cost[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[i][ i_next].T  + cost[s[i]][s[i_next+1]];

        cost_new = seq[0][ i_prev].C                                                    +           //       1st subseq
        seq[i][ i_next].W               * (cost_concat_1) + cost[s[i_next]][s[i]]  +           // concat 2nd subseq
        seq[i_next+1][ dimen].W   * (cost_concat_2) + seq[i_next+1][ dimen].C;   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 seq[0][ i_prev].T       + cost[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1                           + cost[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[i_next][ j_prev].T  + cost[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3                           + cost[s[i]][s[j_next]];

            cost_new = seq[0][ i_prev].C                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            seq[i_next][ j_prev].W      * cost_concat_2 + seq[i_next][ j_prev].C +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            seq[j_next][ dimen].W * cost_concat_4 + seq[j_next][ dimen].C;   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < seq[0][dimen].C - DBL_EPSILON) {
        swap(s, I, J);
        subseq_load(s, seq, I);
        if (cost_best != seq[0][dimen].C) {
            cout << "difere " << endl;
        }

        //return true;
        improve = true;
    }

    //return false;
}

//inline bool search_two_opt(vector<int> & s, tSubseq & seq) {
void search_two_opt(vector<int> & s, tSubseq & seq) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2;
    alignas(DBL_SZ) double cost_best = DBL_MAX;// cost_l1, cost_l2;
    alignas(DBL_SZ) double rev_seq_cost;
    alignas(INT_SZ) int i, j, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = seq[i][ i+1].T;
        for (j = i + 2; j < dimen; ++j) {
            j_next = j + 1;

            rev_seq_cost += cost[s[j-1]][s[j]] * (seq[i][ j].W-1.0);

            cost_concat_1 =                 seq[0][ i_prev].T   + cost[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[i][ j].T        + cost[s[j_next]][s[i]];

            cost_new = seq[0][ i_prev].C                                                        +   //  1st subseq
                    seq[i][ j].W                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
                    seq[j_next][ dimen].W * cost_concat_2 + seq[j_next][ dimen].C;      // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < seq[0][dimen].C - DBL_EPSILON) {
        reverse(s, I, J);
        subseq_load(s, seq);//, I);  // d olho aq hein
        if (cost_best != seq[0][dimen].C) {
            cout << "difere " << endl;
        }

        //return true;
        improve = true;
    }

    //return false;
}

void search_reinsertion(vector<int> & s, tSubseq & seq, const int opt) {
//inline bool search_reinsertion(vector<int> & s, tSubseq & seq, const int opt) {
    alignas(DBL_SZ) double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    alignas(DBL_SZ) double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, k_next, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;
    alignas(INT_SZ) int POS;

    for (i = 1, j = opt +i-1; i < dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

            cost_concat_1 =                 seq[0][ k].T            + cost[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[i][ j].T            + cost[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][ i_prev].T  + cost[s[i_prev]][s[j_next]];

            cost_new = seq[0][ k].C                                                                   +   //       1st subseq
                seq[i][ j].W               * cost_concat_1 + seq[i][ j].C                  +   //  concat 2nd subseq (reinserted seq)
                seq[k_next][ i_prev].W     * cost_concat_2 + seq[k_next][ i_prev].C        +   //  concat 3rd subseq
                seq[j_next][ dimen].W * cost_concat_3 + seq[j_next][ dimen].C;       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < dimen ; ++k) {
            k_next = k + 1;

            cost_concat_1 =                 seq[0][ i_prev].T  + cost[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[j_next][ k].T  + cost[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[i][ j].T       + cost[s[j]][s[k_next]];

            cost_new = seq[0][ i_prev].C                                                                  +   //       1st subseq
                seq[j_next][ k].W          * cost_concat_1 + seq[j_next][ k].C             +   // concat 2nd subseq
                seq[i][ j].W               * cost_concat_2 + seq[i][ j].C                  +   // concat 3rd subseq (reinserted seq)
                seq[k_next][ dimen].W * cost_concat_3 + seq[k_next][ dimen].C;       // concat 4th subseq


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
        subseq_load(s, seq, I < POS+1 ? I : POS+1);

        if (cost_best != seq[0][dimen].C) {
            cout << "difere " << endl;
        }

        improve = true;
        //return true;
    }

    //return false;
}

void RVND(vector<int> & s, tSubseq & seq, tRnd & rnd) {

    vector<string> nome (6);

    nome[SWAP] = "SWAP";
    nome[TWO_OPT] = "TWO_OPT";
    nome[REINSERTION] = "REINSERTION";
    nome[OR_OPT_2] = "OR_OPT_2";
    nome[OR_OPT_3] = "OR_OPT_3";

    alignas(alignof(std::vector<int>)) std::vector<int> neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    alignas(INT_SZ) uint index;
    alignas(INT_SZ) int neighbd;
    //int k = 0;
    bool improve_flag;

    while (!neighbd_list.empty()) {
        //k++;

        index = rand() % neighbd_list.size();
        //cout << info.rnd[info.rnd_index] << endl;
        index = rnd.rnd[rnd.rnd_index++];
        //cout << index << endl;
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";

        //improve_flag = false;
        improve = false;

        switch(neighbd){
            case REINSERTION:
                search_reinsertion(s, seq, REINSERTION);
                //improve_flag = search_reinsertion(s, seq, REINSERTION);
                break;				
            case OR_OPT_2:
                search_reinsertion(s, seq, OR_OPT_2);
                //improve_flag = search_reinsertion(s, seq, OR_OPT_2);
                break;				
            case OR_OPT_3:
                search_reinsertion(s, seq, OR_OPT_3);
                //improve_flag = search_reinsertion(s, seq, OR_OPT_3);
                break;				
            case SWAP:
                search_swap(s, seq);
                //improve_flag = search_swap(s, seq);
                break;
            case TWO_OPT:
                search_two_opt(s, seq);
                //improve_flag = search_two_opt(s, seq);
                break;				
        }

        if (improve) {
            //cout << "improv  " << nome[neighbd] << endl;
            neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
        } else {
            //cout << "delete  " << nome[neighbd] << endl;
            neighbd_list.erase(neighbd_list.begin() + index);
        }

    }
}

std::vector<int> perturb(vector<int> & sl, tRnd & rnd) {
    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    auto s = sl;


    int size_max = std::floor((dimen+1)/10);
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
        A_start = rnd.rnd[rnd.rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        A_end = A_start + rnd.rnd[rnd.rnd_index++];
        //std::cout << "A start  " << A_start << std::endl;
        //std::cout << "A end  " << A_end << std::endl;

        //cout << info.rnd[info.rnd_index] << endl;
        B_start = rnd.rnd[rnd.rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        B_end = B_start + rnd.rnd[rnd.rnd_index++];
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

    return s;
}


void GILS_RVND(int Imax, int Iils, tRnd & rnd) {

    vector<int> solut_partial (dimen+1);
    vector<int> solut_crnt (dimen+1);
    vector<int> solut_best (dimen+1);
    double cost_partial = DBL_MAX;
    double cost_crnt = DBL_MAX;
    double cost_best = DBL_MAX;

    tSubseq seq(dimen+1, vector<seqStruct> (dimen+1));

    for(int i = 0; i < Imax; ++i){
        /**/ int aux = (unsigned)rand() % TABLE_SZ;
        aux = rnd.rnd[rnd.rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	


        solut_crnt = construct(alpha, rnd);
        solut_partial = solut_crnt;

        subseq_load(solut_crnt, seq);
        cost_crnt = seq[0][dimen].C;
        cost_partial = seq[0][dimen].C;

        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", cost_partial);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            RVND(solut_crnt, seq, rnd);
            cost_crnt = seq[0][dimen].C;

            if (cost_crnt < cost_partial - DBL_EPSILON) {
                solut_partial = solut_crnt;
                cost_partial = cost_crnt;
                //solut_partial = solut_crnt;
                iterILS = 0;
            }


            solut_crnt = perturb(solut_partial, rnd);
            subseq_load(solut_crnt, seq);
            //cout << "perturb  " << seq[0][dimen].C << endl;

            iterILS++;
        }

        //subseq_load(solut_partial, info);

        if (cost_partial < cost_best - DBL_EPSILON) {
            solut_best = solut_partial;
            cost_best = cost_partial;
        }

        //after(7);

        //std::cout << "\tCurrent search cost: "<< cost_sl << std::endl;
        std::cout << "\tCurrent best cost: "<< cost_best << std::endl;
        //std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        //std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;
        //std::cout << k << "  Iteracoes " << std::endl;

        std::cout << "SOLUCAO: ";
        for(int i = 0; i < solut_best.size(); i++){
            std::cout << solut_best[i] << " ";
        }
        std::cout << std::endl;

    }
    //std::cout << "Dimension: " << dimension << std::endl;
    printf("COST: %.2lf\n", cost_best);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tRnd rnd;

    std::vector<int> rnd_vec;

    dimen = loadData(&cost, rnd_vec);
    cout << cost[0][dimen-1] << endl;
    rnd.rnd = rnd_vec;
    //print_s(rnd.rnd);
    rnd.rnd_index = 0;


    /*
    for (int i = 0; i <=n; i++) {
        for (int j = i+1; j <=n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */


    srand(clock());

    Iils = dimen < 100 ? dimen : 100;
    auto t1 = high_resolution_clock::now();
    GILS_RVND(Imax, Iils, rnd);
    auto t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();

    double res = (double)duration / 10e2;
    std::cout << "TIME: " << res << std::endl;

    std::cout << "Tamanho RND " << rnd.rnd.size() << std::endl;

    exit(0);
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

