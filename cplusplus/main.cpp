#include <iostream>
#include <cstdint>
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
#define OR_OPT2 	2
#define OR_OPT3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

typedef unsigned uint;

// 3D array
typedef std::vector<std::vector<std::vector<double>>> tSubseq;

typedef struct tInfo {
    double ** cost;
    int dimen;
    uint T;
    uint C;
    uint W;
} tInfo;

typedef struct tSolution {
    std::vector<int> s;
    tSubseq seq;
    double cost;
} tSolution;

tSolution Solution_init(tInfo info) {
    tSolution solut;
    solut.s = vector<int>(info.dimen+1);
    /*
    solut.seq = new double ** [info.dimen+1];
    for (int i = 0; i < info.dimen+1; i++) {
        solut.seq[i] = new double * [info.dimen+1];
        for (int j = 0; j < info.dimen+1; j++) {
            solut.seq[i][j] = new double [3];
        }
    }
    */
    solut.seq = std::vector<std::vector<std::vector<double>>> (
            info.dimen+1, std::vector<std::vector<double>> (
                info.dimen+1, std::vector<double> (3)));
    solut.cost = DBL_MAX;

    return solut;
}

double R_table(int i){
    static const double table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

void print_s(std::vector<int> s) {

    for (int i = 0; i < s.size(); i++)
        std::cout << s[i] << " ";
    std::cout << std::endl;
}

std::vector<int> construct(const double alpha, tInfo info){

    std::vector<int> s = {0};
    std::vector<int> cL(info.dimen-1);
    for(int i = 1; i < info.dimen; ++i){
        cL[i-1] = i;
    }

    int r = 0;
    while (!cL.empty()) {
        std::sort(cL.begin(), cL.end(), 
            [r, &info] (const int i, const int j) {
                return info.cost[i][r] < info.cost[j][r];
            });

        int range = ceil(cL.size() * alpha);
        int index = rand() % range;
        int c = cL[index];
        s.push_back(c);
        r = c;
        print_s(cL);
        cL.erase(cL.begin() + index);
    }

    s.push_back(0);

    return s;
}	

inline void swap_2(std::vector<int> &vec, int i, int j){
    std::iter_swap(vec.begin() + i, vec.begin() + j);
}

inline void two_opt(std::vector<int> &vec, int i, int j){
    std::reverse(vec.begin() + i, vec.begin() + j+1);
}

inline void reinsert(std::vector<int> &vec, int i, int j, int pos){
    if(pos < i){
        std::vector<int> copy;
        copy.reserve(j - i);
        std::vector<int>::iterator t1 = vec.begin() + i;
        std::vector<int>::iterator t2 = vec.begin() + j; 
        copy.insert(copy.begin(),t1 ,t2 );
        vec.erase(t1, t2);
        vec.insert(vec.begin() + pos, copy.begin(), copy.end());
    }else{
        std::vector<int> copy;
        copy.reserve(j - i );
        copy.insert(copy.begin(), vec.begin() + i, vec.begin() + j);
        vec.insert(vec.begin() + pos, copy.begin(), copy.end());
        vec.erase(vec.begin() + i, vec.begin() + j);
    }

}

inline void subseq_load(tSolution & solut, tInfo info){
    alignas(INT_SZ) int i, j, j_prev, k;
    //alignas(INT_SZ) int dim = dimension+1;
    //alignas(1) bool t;
    for (i = 0; i < info.dimen+1; i++) {
        k = 1 - i - (!i);

        solut.seq[i][i][info.T] = 0.0;
        solut.seq[i][i][info.C] = 0.0;
        solut.seq[i][i][info.W] = (double) !(i == 0);
        for (j = i+1; j < info.dimen+1; j++) {
            j_prev = j-1;
            solut.seq[i][j][info.T] = info.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][info.T];
            solut.seq[i][j][info.C] = solut.seq[i][j][info.T] + solut.seq[i][j_prev][info.C];
            solut.seq[i][j][info.W] = j + k;
        }
    }

    solut.cost = solut.seq[0][info.dimen][info.C];
}

/*
inline void subseq_info_load2(std::vector<std::vector<struct subseq>> &seq, std::vector<int> &s, int index = 0){
    alignas(INT_SZ) int i, j, a, k;
    alignas(INT_SZ) int dim = dimension+1;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for(i = 0; i < dim; i++){
        k = 1 - i -(!i);
        t = i == from;
        for(j = from + t; j < dim; j++){
            a = j-1;
            seq[i][j].T = c[s[a]][s[j]] + seq[i][a].T;
            seq[i][j].C = seq[i][j].T + seq[i][a].C;
            seq[i][j].W = j + k ;
        }
        from += t;
    }
}

inline double cost_reverse_calc(std::vector<std::vector<struct subseq>> &seq, std::vector<int> &s, int from, int to){
    double sum = 0;
    for(int i = from; i < to; ++i)
        sum += seq[i][to].T; 
    return sum;
}

inline void neighbor_swap_better(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq){
    alignas(DBL_SZ) double cost, cost1, cost2, cost3, cost4;
    alignas(DBL_SZ) double cost_lower;
    alignas(INT_SZ) int i, j, e, b, d, a, q;
    alignas(INT_SZ) int dim = dimension - 2;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(1) bool t = 1;
    alignas(1) bool l;

    for(i = 1; i < dimension - 1; ++i){
        d = i - 1;
        a = i + 1;
        q = a + 1;

        //consecutive nodes
        cost1 = seq[0][d].T + c[s[d]][s[a]];
        cost2 = cost1 + seq[i][a].T  + c[s[i]][s[q]];

        cost = seq[0][d].C + seq[i][a].W * (cost1) +  c[s[a]][s[i]];
        cost = cost + seq[q][dimension].W * (cost2) +  seq[q][dimension].C;

        if(cost < cost_lower || t){
            cost_lower = cost - DBL_EPSILON;
            i_best = i;
            j_best = a;
            t = 0;
        }

        if(i == dim) continue;

        for(j = q; j < dimension; ++j){
            b = j + 1;
            e = j - 1;

            cost1 = seq[0][d].T + c[s[d]][s[j]];
            cost2 = cost1 + c[s[j]][s[a]];
            cost3 = cost2 + seq[a][e].T + c[s[e]][s[i]];
            cost4 = cost3 + c[s[i]][s[b]];

            cost  = seq[0][d].C + cost1;
            cost += seq[a][e].W * cost2 + seq[a][e].C;
            cost += cost3;
            cost += seq[b][dimension].W * cost4 + seq[b][dimension].C;

            if(cost < cost_lower){
                cost_lower = cost - DBL_EPSILON;
                i_best = i;
                j_best = j;
            }
        }
    }

    if(cost_lower < seq[0][dimension].C - DBL_EPSILON){
        swap_2(s, i_best, j_best);
        subseq_info_load2(seq, s, i_best);
        state = true;
    }
}

inline void neighbor_two_opt_better(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq){
    alignas(DBL_SZ) double cost, cost1, cost2;
    alignas(DBL_SZ) double cost_lower, cost_l1, cost_l2;
    alignas(DBL_SZ) double rev;
    alignas(INT_SZ) int i, j, b, a;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(1) bool t = 1;

    for(i = 1; i < dimension - 1; ++i){
        b = i - 1;

        for(j = i + 2; j < dimension; ++j){
            a = j + 1;
            rev = cost_reverse_calc(seq, s, i, j);

            cost1 = seq[0][b].T + c[s[j]][s[b]];
            cost2 = cost1 + seq[i][j].T  + c[s[a]][s[i]];

            cost_l1 = seq[0][b].C + seq[i][j].W * (cost1) +  rev;
            cost_l2 = seq[a][dimension].W * (cost2) +  seq[a][dimension].C;

            cost = cost_l1 + cost_l2;
            
            if(cost < cost_lower || t){
                cost_lower = cost - DBL_EPSILON;
                i_best = i;
                j_best = j;
                t = 0;

            }
        }
    }

    if(cost_lower < seq[0][dimension].C - DBL_EPSILON){
        two_opt(s, i_best, j_best);
        subseq_info_load2(seq, s, i_best);
        state = true;
    }
}


inline void neighbor_reinsertion_better(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq, int sz){
    alignas(DBL_SZ) double cost, cost1, cost2, cost3;
    alignas(DBL_SZ) double cost_lower, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, g, b, e;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(INT_SZ) int pos_new;
    alignas(INT_SZ) int k_lim = dimension - sz - 1;
    alignas(1) bool l;
    alignas(1) bool t = 1;

    for(i = 1, j = sz + i - 1; i < dimension - sz + 1; ++i, ++j){
        e = j + 1;
        b = i - 1;

        //k -> edges 
        for(k = 0; k < b; ++k){
            g = k + 1;

            cost1 = seq[0][k].T + c[s[k]][s[i]];
            cost2 = cost1 + seq[i][j].T + c[s[j]][s[g]];
            cost3 = cost2 + seq[g][b].T + c[s[b]][s[e]];

            cost_l1 = seq[0][k].C + seq[i][j].W * cost1 + seq[i][j].C; 
            cost_l2 =               seq[g][b].W * cost2 + seq[g][b].C;
            cost_l3 =       seq[e][dimension].W * cost3 + seq[e][dimension].C;

            cost = cost_l1 + cost_l2 + cost_l3;

            if( cost < cost_lower || t){
                cost_lower = cost - DBL_EPSILON;
                i_best = i;
                j_best = j;
                pos_new = k;
                t = 0;
                l = 0;

            }
        }

        for(k = i + sz; k < k_lim; ++k){
            g = k + 1;

            cost1 = seq[0][b].T + c[s[b]][s[e]];
            cost2 =cost1 + seq[e][k].T + c[s[k]][s[i]];
            cost3 = cost2 + seq[i][j].T + c[s[j]][s[g]];

            cost_l1 = seq[0][b].C + seq[e][k].W * cost1 + seq[e][k].C; 
            cost_l2 =               seq[i][j].W * cost2 + seq[i][j].C;
            cost_l3 =       seq[g][dimension].W * cost3 + seq[g][dimension].C;

            cost = cost_l1 + cost_l2 + cost_l3;

            if( cost < cost_lower){
                cost_lower = cost - DBL_EPSILON;
                i_best = i;
                j_best = j;
                pos_new = k;
                l = 1;
            }
        }
    }

    if(cost_lower < seq[0][dimension].C - DBL_EPSILON){
        reinsert(s, i_best, j_best + 1, pos_new + 1);
        int ar[] = {pos_new+1, i_best};
        subseq_info_load2(seq, s, ar[l]);
        state = true;
    }

}

inline void neighbd_list_repopulate(std::vector<int> &list){
    list.clear();
    list = {REINSERTION, OR_OPT_2, OR_OPT_3, SWAP, TWO_OPT};
}

void RVND(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq){

    alignas(alignof(std::vector<int>)) std::vector<int> neighbd_list = {1,2,3,4,5};
    alignas(INT_SZ) int neighbd_rand_index;
    alignas(INT_SZ) int neighbd_rand;
    int k = 0;

    while(!neighbd_list.empty()){
        k++;

        neighbd_rand_index = (unsigned)rand() % neighbd_list.size();
        neighbd_rand = neighbd_list[neighbd_rand_index];

        state = false;

        switch(neighbd_rand){
            case REINSERTION:
                //before();
                neighbor_reinsertion_better(s, seq, REINSERTION);
                //after(REINSERTION);
                break;				
            case OR_OPT2:
                //before();
                neighbor_reinsertion_better(s, seq, OR_OPT2);
                //after(OR_OPT2);
                break;				
            case OR_OPT3:
                //before();
                neighbor_reinsertion_better(s, seq, OR_OPT3);
                //after(OR_OPT3);
                break;				
            case SWAP:
                //before();
                neighbor_swap_better(s, seq);
                //after(SWAP);
                break;
            case TWO_OPT:
                //before();
                neighbor_two_opt_better(s, seq);
                //after(TWO_OPT);
                break;				
        }

        if(state)
            neighbd_list_repopulate(neighbd_list);
        else
            neighbd_list.erase(neighbd_list.begin() + neighbd_rand_index);
        

    }

    //std::cout << k << " RVND iteracoes" << std::endl;
}

void perturb(std::vector<int> &sl, std::vector<int> &s){
}
*/


void GILS_RVND(int Imax, int Iils, tInfo info) {

    tSolution solut_partial = Solution_init(info);
    tSolution solut_crnt = Solution_init(info);
    tSolution solut_best = Solution_init(info);

    for(int i = 0; i < Imax; ++i){
        //before();
        int aux = (unsigned)rand() % TABLE_SZ;
        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	

        solut_crnt.s = construct(alpha, info);
        subseq_load(solut_crnt, info);

        solut_partial = solut_crnt;
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            //RVND(s, subseq_info);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                solut_partial = solut_crnt;
                iterILS = 0;
            }

            //perturb(sl, s);
            //subseq_load(solut, s);
            iterILS++;
        }

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            solut_best = solut_partial;
        }

        //after(7);

        //std::cout << "\tCurrent search cost: "<< cost_sl << std::endl;
        //std::cout << "\tCurrent best cost: "<< cost_final << std::endl;
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
    //printf("COST: %.2lf\n", solut_best.cost);
}

int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tInfo info = {
        .T = 0,
        .C = 1,
        .W = 2
    };


    info.dimen = loadData(&info.cost);

    tSolution solut = Solution_init(info);

    /*
    for (int i = 0; i <=n; i++) {
        for (int j = i+1; j <=n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
    srand(clock());

    solut.s = construct(0.1, info);
    subseq_load(solut, info);

    std::cout << solut.cost << std:: endl;

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

