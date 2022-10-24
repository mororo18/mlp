#ifndef MLP_HPP_
#define MLP_HPP_

#include <vector>
#include <cfloat>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <new>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <vector>

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "tInfo.hpp"
#include "tSolution.hpp"

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define DBL_SZ      8
#define INT_SZ      4
#define TABLE_SZ    26


class MLP {
private:
    tInfo * info;

    std::vector<int> construct(const double alpha, tInfo & info);
    void subseq_load(tSolution & solut, tInfo & info, int index);
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
