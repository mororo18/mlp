#include "MLP.hpp"
#include "Data.hpp"


using std::chrono::steady_clock;


int main(int argc, char **argv){

    tInfo info;


    std::vector<int> rnd;
    double ** cost;
    int dimen = loadData(&cost, rnd);
    info.setDimen(dimen);
    info.setCostPtr(cost);
    info.setRnd(rnd);

    MLP mlp(info);

    mlp.solve();



    return 0;
}

