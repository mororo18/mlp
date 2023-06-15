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
#include "tSolution.hpp"
#include "tInfo.hpp"

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4


/*
typedef struct tSolution {
    std::vector<int> s;
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

    solut.cost = DBL_MAX;

    return solut;
}

static void Solution_cpy( tSolution & src, tSolution & tgt, const tInfo & info) {

    tgt.s = src.s;
    tgt.cost = src.cost;

}

*/

class MLP {

private:

    tInfo info;



    std::vector<int> construct(const double, tInfo &);

    void subseq_load(tSolution & solut, tInfo & info, int index );
    bool search_swap(tSolution & solut, tInfo & info);
    bool search_two_opt(tSolution & solut, tInfo & info);
    bool search_reinsertion(tSolution & solut, tInfo & info, const int opt);
    void RVND(tSolution & solut, tInfo & info);
    std::vector<int> perturb(tSolution * solut, tInfo & info);
    void GILS_RVND(int Imax, int Iils, tInfo & info);

public:
    MLP(tInfo & info);
    void solve();

};

#endif 
