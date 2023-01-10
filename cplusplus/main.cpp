#include "Data.hpp"
#include "MLP.hpp"


using std::chrono::high_resolution_clock;



int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    tInfo info = {};
    info.T = 0;
    info.W = 1;
    info.C = 2;


    std::vector<int> rnd;

    info.dimen = loadData(&info.cost, rnd);
    info.rnd = rnd;
    info.rnd_index = 0;

    MLP mlp(info);

    mlp.solve();


    return 0;
}

