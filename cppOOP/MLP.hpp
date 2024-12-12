#ifndef MLP_HPP
#define MLP_HPP

#include <cstring>
#include <cfloat>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include "tSolution.hpp"
#include "tData.hpp"

#define REINSERTION 1
#define OR_OPT_2 	2
#define OR_OPT_3 	3
#define SWAP 		4
#define TWO_OPT		5
#define TABLE_SZ    26
#define DBL_SZ      8
#define INT_SZ      4

class MLP {

private:

    tData data;

    std::vector<int> construct(const double, tData &);

    void update_subseq_info_matrix(tSolution & solut, tData & data, int index );
    bool search_swap(tSolution & solut, tData & data);
    bool search_two_opt(tSolution & solut, tData & data);
    bool search_reinsertion(tSolution & solut, tData & data, const int opt);
    void RVND(tSolution & solut, tData & data);
    std::vector<int> perturb(tSolution * solut, tData & data);
    void GILS_RVND(int Imax, int Iils, tData & data);

public:
    MLP(tData & data);
    void solve();

};

#endif 
