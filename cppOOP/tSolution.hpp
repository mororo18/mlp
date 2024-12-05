#ifndef _SOLUTION_HPP
#define _SOLUTION_HPP

#include <vector>
#include <cfloat>
#include <algorithm>

#include "tInfo.hpp"

typedef struct tSeqInfo {
    double T;
    double C;
    double W;
} tSeqInfo;


class tSolution {
private:
    std::vector<int> s;
    tSeqInfo ** seq;
    double cost;

public:
    tSolution(tData & info);

    double getCost() {return cost;}
    std::vector<int> getSolutVec() {return s;}

    void setSolutVec(std::vector<int> s) {this->s = s;}

    inline void setT(int i, int j, double T) {this->seq[i][j].T = T;}
    inline void setC(int i, int j, double C) {this->seq[i][j].C = C;}
    inline void setW(int i, int j, double W) {this->seq[i][j].W = W;}

    inline double getT(int i, int j) {return this->seq[i][j].T;}
    inline double getC(int i, int j) {return this->seq[i][j].C;}
    inline double getW(int i, int j) {return this->seq[i][j].W;}

    inline int getPos(int i) {return this->s[i];}
    inline void setCost(double cost) {this->cost = cost;}

    inline void swap(int i, int j){
        std::iter_swap(this->s.begin() + i, this->s.begin() + j);
    }

    inline void reverse(int i, int j){
        std::reverse(this->s.begin() + i, this->s.begin() + j+1);
    }

    inline void reinsert(int i, int j, int pos){
        std::vector<int> seq (this->s.begin() + i, this->s.begin() +j+1);
        if(pos < i){
            this->s.erase(this->s.begin() + i, this->s.begin() + j+1);
            this->s.insert(this->s.begin() + pos, seq.begin(), seq.end());
        }else{
            this->s.insert(this->s.begin() + pos, seq.begin(), seq.end());
            this->s.erase(this->s.begin() + i, this->s.begin() + j+1);
        }
    }


    void copy(tSolution &);
    void print();
};
#endif
