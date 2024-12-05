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

#define MATRIX

using std::chrono::high_resolution_clock;

typedef unsigned uint;

// 3D array
typedef std::vector<std::vector<std::vector<double>>> tSubseq;

typedef struct tData {
    double ** cost;
    int dimen;
    uint T;
    uint C;
    uint W;
    vector<int> rnd;
    uint rnd_index;
} tData;

typedef struct tSolution {
    std::vector<int> s;
    //tSeqInfo ** seq;
    double *** seq;
    double cost;
} tSolution;

tSolution Solution_init(tData data) {
    tSolution solut;
    solut.s = vector<int>(data.dimen+1);

    solut.seq = new double ** [data.dimen+1];
    for (int i = 0; i < data.dimen+1; i++) {
        solut.seq[i] = new double * [data.dimen+1];
        for (int j = 0; j < data.dimen+1; j++) {
            solut.seq[i][j] = new double [3];
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

inline void Solution_cpy( tSolution & src, tSolution & tgt, const tData & data) {

    tgt.s = src.s;
    tgt.cost = src.cost;

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

void print_s(std::vector<int> s) {

    for (int i = 0; i < s.size(); i++)
        std::cout << s[i]+1 << " ";
    std::cout << std::endl;
}

int partition(std::vector<int>& arr, int left, int right, const tData& data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (data.cost[r][arr[j]] < data.cost[r][pivot]) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

void quicksort(std::vector<int>& arr, int left, int right, const tData& data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

void sort(std::vector<int>& arr, int r, const tData& data) {
    quicksort(arr, 0, arr.size() - 1, data, r);
}

std::vector<int> construct(const double alpha, tData & data){

    std::vector<int> s = {0};
    s.reserve(data.dimen+1);
    std::vector<int> cL(data.dimen-1);
    for(int i = 1; i < data.dimen; ++i){
        cL[i-1] = i;
    }

    int r = 0;
    while (!cL.empty()) {
        sort(cL, r, data);

        /**/
        int range = std::ceil(cL.size() * alpha);
        int index = range > 0 ? rand() % range : 0;
        /**/

        //std::cout << data.rnd[data.rnd_index]<< std::endl;
        index = data.rnd[data.rnd_index++];
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



inline void update_subseq_info_matrix(tSolution & solut, const tData & data, int index = 0){
    alignas(INT_SZ) int i, j, j_prev, k;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for (i = 0; i < data.dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut.seq[i][i][data.T] = 0.0;
        solut.seq[i][i][data.C] = 0.0;
        solut.seq[i][i][data.W] = (double) !(i == 0);

        for (j = i+1; j < data.dimen+1; j++) {
            j_prev = j-1;
            
            solut.seq[i][j][data.T] = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][data.T];
            solut.seq[i][j][data.C] = solut.seq[i][j][data.T] + solut.seq[i][j_prev][data.C];
            solut.seq[i][j][data.W] = j + k;

        }
        from += t;
    }

    solut.cost = solut.seq[0][data.dimen][data.C];
}

inline bool search_swap(tSolution & solut, const tData & data) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    alignas(DBL_SZ) double cost_best = DBL_MAX;
    alignas(INT_SZ) int i, j, j_prev, j_next, i_prev, i_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data.dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut.seq[0][i_prev][data.T] + data.cost[solut.s[i_prev]][solut.s[i_next]];
        cost_concat_2 = cost_concat_1 + solut.seq[i][i_next][data.T]  + data.cost[solut.s[i]][solut.s[i_next+1]];

        cost_new = solut.seq[0][i_prev][data.C]                                                    +           //       1st subseq
        solut.seq[i][i_next][data.W]               * (cost_concat_1) + data.cost[solut.s[i_next]][solut.s[i]]  +           // concat 2nd subseq
        solut.seq[i_next+1][data.dimen][data.W]   * (cost_concat_2) + solut.seq[i_next+1][data.dimen][data.C];   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data.dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut.seq[0][i_prev][data.T]       + data.cost[solut.s[i_prev]][solut.s[j]];
            cost_concat_2 = cost_concat_1                           + data.cost[solut.s[j]][solut.s[i_next]];
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][j_prev][data.T]  + data.cost[solut.s[j_prev]][solut.s[i]];
            cost_concat_4 = cost_concat_3                           + data.cost[solut.s[i]][solut.s[j_next]];

            cost_new = solut.seq[0][i_prev][data.C]                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut.seq[i_next][j_prev][data.W]      * cost_concat_2 + solut.seq[i_next][j_prev][data.C] +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut.seq[j_next][data.dimen][data.W] * cost_concat_4 + solut.seq[j_next][data.dimen][data.C];   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        swap(solut.s, I, J);
        update_subseq_info_matrix(solut, data, I);
        return true;
    }

    return false;
}

inline bool search_two_opt(tSolution & solut, const tData & data) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2;
    alignas(DBL_SZ) double cost_best = DBL_MAX;// cost_l1, cost_l2;
    alignas(DBL_SZ) double rev_seq_cost;
    alignas(INT_SZ) int i, j, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data.dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut.seq[i][i+1][data.T];
        for (j = i + 2; j < data.dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += data.cost[solut.s[j-1]][solut.s[j]] * (solut.seq[i][j][data.W]-1.0);

          cost_concat_1 =                 solut.seq[0][i_prev][data.T]   + data.cost[solut.s[j]][solut.s[i_prev]];
          cost_concat_2 = cost_concat_1 + solut.seq[i][j][data.T]        + data.cost[solut.s[j_next]][solut.s[i]];

          cost_new = solut.seq[0][i_prev][data.C]                                                        +   //  1st subseq
              solut.seq[i][j][data.W]                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut.seq[j_next][data.dimen][data.W] * cost_concat_2 + solut.seq[j_next][data.dimen][data.C];      // concat 3rd subseq

            
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
        return true;
    }

    return false;
}

inline bool search_reinsertion(tSolution & solut, const tData & data, const int opt) {
    alignas(DBL_SZ) double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    alignas(DBL_SZ) double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, k_next, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;
    alignas(INT_SZ) int POS;

    for (i = 1, j = opt +i-1; i < data.dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 solut.seq[0][k][data.T]            + data.cost[solut.s[k]][solut.s[i]];
          cost_concat_2 = cost_concat_1 + solut.seq[i][j][data.T]            + data.cost[solut.s[j]][solut.s[k_next]];
          cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev][data.T]  + data.cost[solut.s[i_prev]][solut.s[j_next]];

          cost_new = solut.seq[0][k][data.C]                                                                   +   //       1st subseq
              solut.seq[i][j][data.W]               * cost_concat_1 + solut.seq[i][j][data.C]                  +   //  concat 2nd subseq (reinserted seq)
              solut.seq[k_next][i_prev][data.W]     * cost_concat_2 + solut.seq[k_next][i_prev][data.C]        +   //  concat 3rd subseq
              solut.seq[j_next][data.dimen][data.W] * cost_concat_3 + solut.seq[j_next][data.dimen][data.C];       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data.dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut.seq[0][i_prev][data.T]  + data.cost[solut.s[i_prev]][solut.s[j_next]];
          cost_concat_2 = cost_concat_1 + solut.seq[j_next][k][data.T]  + data.cost[solut.s[k]][solut.s[i]];
          cost_concat_3 = cost_concat_2 + solut.seq[i][j][data.T]       + data.cost[solut.s[j]][solut.s[k_next]];

          cost_new = solut.seq[0][i_prev][data.C]                                                                  +   //       1st subseq
                  solut.seq[j_next][k][data.W]          * cost_concat_1 + solut.seq[j_next][k][data.C]             +   // concat 2nd subseq
                  solut.seq[i][j][data.W]               * cost_concat_2 + solut.seq[i][j][data.C]                  +   // concat 3rd subseq (reinserted seq)
                  solut.seq[k_next][data.dimen][data.W] * cost_concat_3 + solut.seq[k_next][data.dimen][data.C];       // concat 4th subseq
          
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
        update_subseq_info_matrix(solut, data, I < POS+1 ? I : POS+1);
        return true;
    }

    return false;
}


void RVND(tSolution & solut, tData & data) {

    alignas(alignof(std::vector<int>)) std::vector<int> neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    alignas(INT_SZ) uint index;
    alignas(INT_SZ) int neighbd;
    //int k = 0;
    bool improve_flag;

    //cout << "RVND" << endl;
    while (!neighbd_list.empty()) {
        //k++;

        index = rand() % neighbd_list.size();
        //cout << data.rnd[data.rnd_index] << endl;
        index = data.rnd[data.rnd_index++];
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";

        improve_flag = false;

        switch(neighbd){
            case REINSERTION:
                //before();
                improve_flag = search_reinsertion(solut, data, REINSERTION);
                //after(REINSERTION);
                break;				
            case OR_OPT_2:
                //before();
                improve_flag = search_reinsertion(solut, data, OR_OPT_2);
                //after(OR_OPT2);
                break;				
            case OR_OPT_3:
                //before();
                improve_flag = search_reinsertion(solut, data, OR_OPT_3);
                //after(OR_OPT3);
                break;				
            case SWAP:
                //before();
                improve_flag = search_swap(solut, data);
                //after(SWAP);
                break;
            case TWO_OPT:
                //before();
                improve_flag = search_two_opt(solut, data);
                //after(TWO_OPT);
                break;				
        }
        //std::cout << (improve_flag ? "True" : "False") << std::endl;
        if (improve_flag) {
            neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
        } else {
            //std::cout << index << "  " << neighbd_list.size() << std::endl;
            //std::cout << solut.cost << std::endl;
            
            //std::cout << data.rnd_index << std::endl;
            neighbd_list.erase(neighbd_list.begin() + index);
        }

        //std::cout << "cost  " << solut.cost << std::endl ;


    }

    //exit(0);
    //std::cout << k << " RVND iteracoes" << std::endl;
}

std::vector<int> perturb(tSolution * solut, tData & data) {
    auto s = solut->s;
    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = std::floor((data.dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;
    //std::cout << "perturbing\n";
    //print_s(s);
    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        /**/
        int max = (data.dimen+1) -2 -size_max;
        A_start = rand() % max + 1;
        A_end = A_start + rand() % (size_max - size_min + 1) + size_min;

        B_start = rand() % max + 1;
        B_end = B_start + rand() % (size_max - size_min + 1) + size_min;
        /**/



        //std::cout << "paa\n";

        //cout << data.rnd[data.rnd_index] << endl;
        A_start = data.rnd[data.rnd_index++];
        //cout << data.rnd[data.rnd_index] << endl;
        A_end = A_start + data.rnd[data.rnd_index++];
        //std::cout << "A start  " << A_start << std::endl;
        //std::cout << "A end  " << A_end << std::endl;

        //cout << data.rnd[data.rnd_index] << endl;
        B_start = data.rnd[data.rnd_index++];
        //cout << data.rnd[data.rnd_index] << endl;
        B_end = B_start + data.rnd[data.rnd_index++];
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
    //update_subseq_info_matrix(solut, data);

    return s;
}


void GILS_RVND(int Imax, int Iils, tData & data) {

    tSolution solut_partial = Solution_init(data);
    tSolution solut_crnt = Solution_init(data);
    tSolution solut_best = Solution_init(data);

    for(int i = 0; i < Imax; ++i){
        /**/ int aux = (unsigned)rand() % TABLE_SZ;
        aux = data.rnd[data.rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	


        solut_crnt.s = construct(alpha, data);
        //print_s(solut_crnt.s);
        update_subseq_info_matrix(solut_crnt, data);

        //solut_partial = solut_crnt;
        Solution_cpy(solut_crnt, solut_partial, data);
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            RVND(solut_crnt, data);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                Solution_cpy(solut_crnt, solut_partial, data);
                //solut_partial = solut_crnt;
                iterILS = 0;
            }

            solut_crnt.s = perturb(&solut_partial, data);
            update_subseq_info_matrix(solut_crnt, data);
            //exit(0);
            //std::cout << "ITER  " << iterILS << std::endl;
            iterILS++;
        }

        //update_subseq_info_matrix(solut_partial, data);

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(solut_partial, solut_best, data);
            //solut_best = solut_partial;
        }

        //after(7);

        //std::cout << "\tCurrent search cost: "<< cost_sl << std::endl;
        std::cout << "\tCurrent best cost: "<< solut_best.cost << std::endl;
        //std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        //std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;
        //std::cout << k << "  Iteracoes " << std::endl;

        std::cout << "SOLUCAO: ";
        for(int i = 0; i < solut_best.s.size(); i++){
            std::cout << solut_best.s[i] << " ";
        }
        std::cout << std::endl;

    }
    //std::cout << "Dimension: " << dimension << std::endl;
    printf("COST: %.2lf\n", solut_best.cost);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tData data = {};
    data.T = 0;
    data.W = 1;
    data.C = 2;


    std::vector<int> rnd;

    data.dimen = loadData(&data.cost, rnd);
    data.rnd = rnd;
    //print_s(rnd);
    data.rnd_index = 0;

    tSolution solut = Solution_init(data);

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
    auto t1 = high_resolution_clock::now();
    GILS_RVND(Imax, Iils, data);
    auto t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();

    double res = (double)duration / 10e2;
    std::cout << "TIME: " << res << std::endl;

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

