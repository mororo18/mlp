#include "readData.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char ** argv){

    int dimension;
    double ** cost;
    
    readData(argc, argv, &dimension, &cost);

    ofstream file("distance_matrix");

    file << dimension << endl;
    for(int i = 1; i <= dimension; i++){
        for(int j = i+1; j <= dimension; j++){
            file << cost[i][j] << " ";
            if (cost[i][j] == 0.0) {
                cout << "yess\n";
                cout << i << " " << j << endl;
            }
        }
        if(i != dimension)
            file << endl;
        else
            file << "EOF" << endl;;

    }

    string inst(argv[1]+10);

    file << inst << std::endl ;
    //cout << inst.substr(0, inst.find('.')) << endl;
    //cout <<    << endl;

    file << "RND\n";

    fstream r_values("rand_iter_values/" + inst.substr(0, inst.find('.')) + ".rnd", ios_base::in);

    //cout << "Good state: " << r_values.good() << endl;
    if (r_values.good() == false) {
        cout << "Aborted: Ainda nao arquivo \'.rnd\' para essa instancia.\n"; 
        exit(0);
    }

    string line;

    int n = 0;
    while (getline(r_values, line)) {
        file << line << endl;
    }


    file.close();

    return 0;
}
