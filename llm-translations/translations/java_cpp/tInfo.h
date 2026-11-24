#ifndef TINFO_H
#define TINFO_H

#include <vector>

class tInfo {
private:
    std::vector<int>    rnd;
    int                 rnd_index;

    int                 dimension;
    std::vector<std::vector<double>> c;

public:
    static constexpr double EPSILON = 1e-16;

    static constexpr int T = 0;
    static constexpr int C = 1;
    static constexpr int W = 2;

    static constexpr int SWAP        = 0;
    static constexpr int REINSERTION = 1;
    static constexpr int OR_OPT2     = 2;
    static constexpr int OR_OPT3     = 3;
    static constexpr int TWO_OPT     = 4;

    tInfo(int dimen, const std::vector<std::vector<double>> & cost, const std::vector<int> & rnd_arr) {
        rnd = rnd_arr;
        dimension = dimen;
        c = cost;
        rnd_index = 0;
    }

    int getDimen()                  {return dimension;}
    double getCost(int i, int j)    {return c[i][j];}
    int rndCrnt()                   {return rnd[rnd_index++];}
};

#endif
