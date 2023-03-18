#include "MLP.hpp"
#include "Data.hpp"


using std::chrono::steady_clock;


int main(int argc, char **argv){

    Data data;
    data.load();

    tInfo info;
    info.setDimen(data.getDimen());
    info.setCostPtr(data.getCostPtr());
    info.setRnd(data.getRndVec());

    MLP mlp(info);

    mlp.solve();

    return 0;
}

