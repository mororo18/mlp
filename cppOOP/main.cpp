#include <cstring>

#include "Data.hpp"
#include "MLP.hpp"

int main(int argc, char **argv){
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        }
    }

    std::vector<int> rnd;
    double ** cost;
    int dimen;

    dimen = loadData(&cost, rnd);

    tData data;
    data.setDimen(dimen);
    data.setCostPtr(cost);
    data.setRnd(rnd);

    MLP mlp(data);

    mlp.solve(verbose);

    return 0;
}

