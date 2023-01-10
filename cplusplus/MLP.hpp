#ifndef MLP_HPP
#define MLP_HPP

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

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

typedef struct tInfo {
    double ** cost;
    int dimen;
    int T;
    int C;
    int W;
    std::vector<int> rnd;
    int rnd_index;
} tInfo;

typedef struct tSolution {
    std::vector<int> s;
    //tSeqInfo ** seq;
    double *** seq;
    double cost;
} tSolution;

static
tSolution Solution_init(tInfo info) {
    tSolution solut;
    solut.s = std::vector<int>(info.dimen+1);

    solut.seq = new double ** [info.dimen+1];
    for (int i = 0; i < info.dimen+1; i++) {
        solut.seq[i] = new double * [info.dimen+1];
        for (int j = 0; j < info.dimen+1; j++) {
            solut.seq[i][j] = new double [3];
        }
    }

    /*
    solut.seq = new tSeqInfo * [info.dimen+1];
    for (int i = 0; i < info.dimen+1; i++) {
        solut.seq[i] = new tSeqInfo [info.dimen+1];
        memset(solut.seq[i], 0.0, (info.dimen+1)*sizeof(tSeqInfo));
    }
    */

    solut.cost = DBL_MAX;

    return solut;
}

static void Solution_cpy( tSolution & src, tSolution & tgt, const tInfo & info) {

    tgt.s = src.s;
    tgt.cost = src.cost;

    /*
    for (int i = 0; i < info.dimen+1; i++) {
        for (int j = 0; j < info.dimen+1; j++) {
            //memcpy(tgt.seq[i][j], src.seq[i][j], 3 * sizeof(double));
            std::copy(src.seq[i][j], src.seq[i][j] + 3, tgt.seq[i][j]);
        }
    }
    */

}

class MLP {

private:

    tInfo info;



    std::vector<int> construct(const double, tInfo &);

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

    void subseq_load(tSolution & solut, const tInfo & info, int index );
    bool search_swap(tSolution & solut, const tInfo & info);
    bool search_two_opt(tSolution & solut, const tInfo & info);
    bool search_reinsertion(tSolution & solut, const tInfo & info, const int opt);
    void RVND(tSolution & solut, tInfo & info);
    std::vector<int> perturb(tSolution * solut, tInfo & info);
    void GILS_RVND(int Imax, int Iils, tInfo & info);

public:
    MLP(tInfo & info);
    void solve();

};

#endif 
