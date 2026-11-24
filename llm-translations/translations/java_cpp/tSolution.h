#ifndef TSOLUTION_H
#define TSOLUTION_H

#include "tArray.h"

#include <vector>

class tSolution {
public:
    static constexpr bool TEST = true;

private:
    tArray s;

    std::vector<std::vector<std::vector<double>>> seq3d;
    std::vector<std::vector<double>> seq2d;

    double cost;
    static constexpr int Z = 3;

public:
    tSolution(int dimen, double c) {
if (TEST)
        seq2d.assign(dimen+1, std::vector<double>((dimen+1)*Z, 0.0));
else
        seq3d.assign(dimen+1, std::vector<std::vector<double>>(dimen+1, std::vector<double>(Z, 0.0)));
        cost = c;
    }

    void storeSolut(tArray sN) {s = std::move(sN);}

    double getSeq(int i, int j, int k) {
if (TEST)
        return seq2d[i][j*Z + k];
else
        return seq3d[i][j][k];
    }

    void setSeq(int i, int j, int k, double value) {
if (TEST)
        seq2d[i][j*Z + k] = value;
else
        seq3d[i][j][k] = value;
    }

    double getCost() {return cost;}
    void setCost(double value) {cost = value;}

    tArray & getSolut() {return s;}
    tArray getSolutCpy() {return s.getArrayCopy();}
};

#endif
