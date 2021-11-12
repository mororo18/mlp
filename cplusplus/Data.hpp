#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double ** matrix_allocate(size_t sz) {
    double ** ptr = new double*[sz];

    for ( int i = 0; i < sz; i++ ) {
        ptr[i] = new double [sz];
    }
    return ptr;
}

int loadData(double *** matrix) {
    int dimension;
    fstream file;

    file.open("../distance_matrix", ios::in);

    if (file.is_open()) {
        string line;

        getline(file, line);
        dimension = stoi(line);
        *matrix = matrix_allocate((size_t) dimension);
        //cout << "dimension " << dimension << endl;

        int i = 0;
        while (getline(file, line)) {
            //cout << line << "\n";
            int j = i + 1;
            while (line.find(" ") != string::npos) {
                size_t pos = line.find(" ");
                double value = stod(line.substr(0, pos));
                line = line.substr(pos+1);
                (*matrix)[i][j] = value;
                (*matrix)[j][i] = value;
                //cout << line << endl;
                j++;
            }

            i++;
        }

    }

    file.close();
    return dimension;
}
