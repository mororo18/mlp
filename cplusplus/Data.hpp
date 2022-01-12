#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

double ** matrix_allocate(size_t sz) {
    double ** ptr = new double*[sz];

    for ( int i = 0; i < sz; i++ ) {
        ptr[i] = new double [sz];
    }
    return ptr;
}

int loadData(double *** matrix, vector<int> & rnd) {
    int dimension;
    fstream file;

    file.open("../distance_matrix", ios::in);

    if (file.is_open()) {
        string line;

        getline(file, line);
        dimension = stoi(line);
        *matrix = matrix_allocate((size_t) dimension);
        //cout << "dimension " << dimension << endl;

        bool flag = false;
        int rnd_size = -1;
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


            if (flag) {

                if (rnd_size < 0) {
                    rnd_size = stoi(line);
                    rnd.reserve(rnd_size);
                } else {
                    cout << line << endl;
                    rnd.push_back(stoi(line));
                }
            }

            if (line == "RND") flag = true;
        }

    }

    file.close();
    return dimension;
}
