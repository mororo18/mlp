#ifndef _INFO_HPP
#define _INFO_HPP

#include <vector>

class tInfo {
private:
    double ** cost;
    int dimension;
    std::vector<int> rnd;
    int rnd_index;
public:
    tInfo();    
    inline void setRnd(std::vector<int> rnd)   {this->rnd          = rnd; 
                                        this->rnd_index     = 0;}
    inline void setDimen(int dimen)            {this->dimension    = dimen;}
    inline void setCostPtr(double ** ptr)      {this->cost         = ptr;}


    inline int getRndCrnt()                    {return this->rnd[this->rnd_index++];}
    inline int getDimen()                      {return this->dimension;}
    inline int getCost(int i, int j)           {return this->cost[i][j];}


    static const int T = 0;
    static const int C = 1;
    static const int W = 2;

};
#endif
