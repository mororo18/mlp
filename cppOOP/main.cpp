#include "MLP.hpp"
#include "Data.hpp"

int main(int argc, char **argv){

    Data instance;
    instance.load();

    tData data;
    data.setDimen(instance.getDimen());
    data.setCostPtr(instance.getCostPtr());
    data.setRnd(instance.getRndVec());

    MLP mlp(data);

    mlp.solve();

    return 0;
}

