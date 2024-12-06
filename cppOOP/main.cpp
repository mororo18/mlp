#include "Data.hpp"
#include "MLP.hpp"

int main(int argc, char **argv){
    std::vector<int> rnd;
    double ** cost;
    int dimen;

    dimen = loadData(&cost, rnd);

    tData data;
    data.setDimen(dimen);
    data.setCostPtr(cost);
    data.setRnd(rnd);

    MLP mlp(data);

    mlp.solve();

    return 0;
}

