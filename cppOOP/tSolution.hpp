#ifndef _SOLUTION_HPP
#define _SOLUTION_HPP

#include <vector>
#include <cfloat>
#include <algorithm>
#include "tData.hpp"

enum {INDEX_T, INDEX_C, INDEX_W};

class tSolution {
private:
    std::vector<int> s;
    double *** seq;
    double cost;
    tData * data;

    //void setCost(double c) {this->cost = c;}
public:
    tSolution(tData & data);

    double getCost() {return cost;}
    std::vector<int> getSolutVec() {return s;}

    void setSolutVec(std::vector<int> s) {this->s = s;}

    inline void setT(int i, int j, double T) {this->seq[i][j][INDEX_T] = T;}
    inline void setC(int i, int j, double C) {this->seq[i][j][INDEX_C] = C;}
    inline void setW(int i, int j, double W) {this->seq[i][j][INDEX_W] = W;}

    inline double getT(int i, int j) {return this->seq[i][j][INDEX_T];}
    inline double getC(int i, int j) {return this->seq[i][j][INDEX_C];}
    inline double getW(int i, int j) {return this->seq[i][j][INDEX_W];}

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

    void validate();
    double recalcCost();
};
#endif
