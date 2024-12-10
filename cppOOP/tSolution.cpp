#include "tSolution.hpp"

#include <iostream>

tSolution::tSolution(tData & data) : s(data.getDimen()+1) {
    this->seq = new double ** [data.getDimen()+1];
    for (int i = 0; i < data.getDimen()+1; i++) {
        this->seq[i] = new double * [data.getDimen()+1];
        for (int j = 0; j < data.getDimen()+1; j++) {
            this->seq[i][j] = new double [3];
        }
    }

    this->cost = DBL_MAX;

    this->data = &data;
}

void tSolution::copy(tSolution & solut) {
    this->s = solut.getSolutVec();
    this->cost = solut.getCost();
}

void tSolution::print() {
    for (int i = 0; i < this->s.size(); i++)
        std::cout << this->s[i] << " ";
    std::cout << std::endl;
}

void tSolution::validate() {	
    std::vector<bool>  is (this->s.size()-1);
    for (int i = 0; i < this->s.size(); i++) {
        is[this->s[i]] = true;
    }

    for (int i = 0; i < is.size(); i++) {
        if (is[i] != true) {
            std::cout << "Invalido\n";
            exit(0);
        }
    }
}

double tSolution::recalcCost() {
    double total = 0;
    int n = this->s.size()-1;
    for (int i = 0; i < this->s.size()-1; i++) {
        total += data->getCost(getPos(i), getPos(i+1)) * n--;
    }
    return total;
}
