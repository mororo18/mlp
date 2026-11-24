#ifndef DATA_H
#define DATA_H

#include <vector>

class Data {
private:
    int dimension;
    std::vector<std::vector<double>> matrix;
    std::vector<int> rnd;
    int rnd_size;

public:
    Data();

    int                     getDimension()              {return dimension;}
    double                  getDistance(int i, int j)   {return matrix[i][j];}
    std::vector<int> &      getRnd()                    {return rnd;}
    int                     getRndSize()                {return rnd_size;}

    void loadData();
};

#endif
