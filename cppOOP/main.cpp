#include <cstring>

#include "MLP.hpp"
#include "Data.hpp"

int main(int argc, char **argv){

    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    Data instance;
    instance.load();

    tData data;
    data.setDimen(instance.getDimen());
    data.setCostPtr(instance.getCostPtr());
    data.setRnd(instance.getRndVec());

    MLP mlp(data);

    mlp.solve(verbose);

    return 0;
}

