#include "Data.hpp"
#include "MLP.hpp"


using std::chrono::high_resolution_clock;



int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    std::vector<int> rnd;
    double ** cost;
    int dimen;

    dimen = loadData(&cost, rnd);

    tInfo info;
    info.setDimen(dimen);
    info.setCostPtr(cost);
    info.setRnd(rnd);

    MLP mlp(info);

    mlp.solve();


    return 0;
}

