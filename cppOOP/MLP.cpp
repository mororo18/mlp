#include "MLP.hpp"

#include <ctime>
#include <iostream>

MLP::MLP(tData & data) {
    this->data = data;
}

void MLP::solve() {
    tSolution solut(data);

    int Imax = 10;
    int Iils;
    Iils = data.getDimen() < 100 ? data.getDimen() : 100;

    size_t start = clock();
    GILS_RVND(Imax, Iils, data);
    double cpu_time = (double)(clock() - start) / CLOCKS_PER_SEC ;

    std::cout << "TIME: " << cpu_time << std::endl;
}


static
double R_table(int i){
    static const double table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    return table[i];
}

static
void print_s(std::vector<int> s) {

    for (int i = 0; i < s.size(); i++)
        std::cout << s[i]+1 << " ";
    std::cout << std::endl;
}

static
int partition(std::vector<int>& arr, int left, int right, tData& data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (data.getCost(r, arr[j]) < data.getCost(r, pivot)) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

static
void quicksort(std::vector<int>& arr, int left, int right, tData& data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

static
void sort(std::vector<int>& arr, int r, tData& data) {
    quicksort(arr, 0, arr.size() - 1, data, r);
}

std::vector<int> MLP::construct(const double alpha, tData & data){

    std::vector<int> s = {0};
    s.reserve(data.getDimen()+1);
    std::vector<int> cL(data.getDimen()-1);
    for(int i = 1; i < data.getDimen(); ++i){
        cL[i-1] = i;
    }

    int r = 0;
    while (!cL.empty()) {
        sort(cL, r, data);

        int index = data.getRndCrnt();
        int c = cL[index];
        s.push_back(c);
        r = c;
        cL.erase(cL.begin() + index);
    }

    s.push_back(0);

    return s;
}	

void MLP::update_subseq_info_matrix(tSolution & solut, tData & data, int index = 0){
    alignas(INT_SZ) int i, j, j_prev, k;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for (i = 0; i < data.getDimen()+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut.setT(i, i, 0.0);
        solut.setC(i, i, 0.0);
        solut.setW(i, i,(double) !(i == 0));

        for (j = i+1; j < data.getDimen()+1; j++) {
            j_prev = j-1;
            
            double T = data.getCost(solut.getPos(j_prev), solut.getPos(j)) + solut.getT(i, j_prev);
            solut.setT(i, j, T); 

            double C = solut.getT(i, j) + solut.getC(i, j_prev);
            solut.setC(i, j, C);

            double W = j + k;
            solut.setW(i, j, W);

        }
        from += t;
    }

    solut.setCost(solut.getC(0, data.getDimen()));
}

bool MLP::search_swap(tSolution & solut, tData & data) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    alignas(DBL_SZ) double cost_best = DBL_MAX;
    alignas(INT_SZ) int i, j, j_prev, j_next, i_prev, i_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data.getDimen()-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

	
        //consecutive nodes
        cost_concat_1 =                 solut.getT(0, i_prev) + data.getCost(solut.getPos(i_prev), solut.getPos(i_next));
        cost_concat_2 = cost_concat_1 + solut.getT(i, i_next)  + data.getCost(solut.getPos(i), solut.getPos(i_next+1));

        cost_new = solut.getC(0, i_prev)                                                    +           //       1st subseq
        solut.getW(i, i_next)               * (cost_concat_1) + data.getCost(solut.getPos(i_next), solut.getPos(i))  +           // concat 2nd subseq
        solut.getW(i_next+1, data.getDimen())   * (cost_concat_2) + solut.getC(i_next+1, data.getDimen());   // concat 3rd subseq


        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < data.getDimen(); ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut.getT(0, i_prev)       + data.getCost(solut.getPos(i_prev), solut.getPos(j));
            cost_concat_2 = cost_concat_1                           + data.getCost(solut.getPos(j), solut.getPos(i_next));
            cost_concat_3 = cost_concat_2 + solut.getT(i_next, j_prev)  + data.getCost(solut.getPos(j_prev), solut.getPos(i));
            cost_concat_4 = cost_concat_3                           + data.getCost(solut.getPos(i), solut.getPos(j_next));

            cost_new = solut.getC(0, i_prev)                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut.getW(i_next, j_prev)      * cost_concat_2 + solut.getC(i_next, j_prev) +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut.getW(j_next, data.getDimen()) * cost_concat_4 + solut.getC(j_next, data.getDimen());   // concat 5th subseq


            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.getCost() - DBL_EPSILON) {
        solut.swap(I, J);
        update_subseq_info_matrix(solut, data, I);

        return true;
    }

    return false;
}

bool MLP::search_two_opt(tSolution & solut, tData & data) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2;
    alignas(DBL_SZ) double cost_best = DBL_MAX;// cost_l1, cost_l2;
    alignas(DBL_SZ) double rev_seq_cost;
    alignas(INT_SZ) int i, j, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < data.getDimen()-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut.getT(i, i+1);
        for (j = i + 2; j < data.getDimen(); ++j) {
            j_next = j + 1;



          rev_seq_cost += data.getCost(solut.getPos(j-1), solut.getPos(j)) * (solut.getW(i, j)-1.0);


          cost_concat_1 =                 solut.getT(0, i_prev)   + data.getCost(solut.getPos(j), solut.getPos(i_prev));
          cost_concat_2 = cost_concat_1 + solut.getT(i, j)        + data.getCost(solut.getPos(j_next), solut.getPos(i));

          cost_new = solut.getC(0, i_prev)                                                        +   //  1st subseq
              solut.getW(i, j)                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut.getW(j_next, data.getDimen()) * cost_concat_2 + solut.getC(j_next, data.getDimen());      // concat 3rd subseq

   
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.getCost() - DBL_EPSILON) {
        solut.reverse(I, J);
        update_subseq_info_matrix(solut, data);

        solut.validate();

        if (solut.getCost() != cost_best) {
            std::cout << "deu merda 2\n";

            std::cout << "Cost best = " << cost_best << "\nCost calc = " << solut.getCost() << "\n";

        }

        return true;
    }

    return false;
}

bool MLP::search_reinsertion(tSolution & solut, tData & data, const int opt) {
    alignas(DBL_SZ) double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    alignas(DBL_SZ) double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, k_next, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;
    alignas(INT_SZ) int POS;

    for (i = 1, j = opt +i-1; i < data.getDimen()-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;


          cost_concat_1 =                 solut.getT(0, k)            + data.getCost(solut.getPos(k), solut.getPos(i));
          cost_concat_2 = cost_concat_1 + solut.getT(i, j)            + data.getCost(solut.getPos(j), solut.getPos(k_next));
          cost_concat_3 = cost_concat_2 + solut.getT(k_next, i_prev)  + data.getCost(solut.getPos(i_prev), solut.getPos(j_next));

          cost_new = solut.getC(0, k)                                                                   +   //       1st subseq
              solut.getW(i, j)               * cost_concat_1 + solut.getC(i, j)                  +   //  concat 2nd subseq (reinserted seq)
              solut.getW(k_next, i_prev)     * cost_concat_2 + solut.getC(k_next, i_prev)        +   //  concat 3rd subseq
              solut.getW(j_next, data.getDimen()) * cost_concat_3 + solut.getC(j_next, data.getDimen());       // concat 4th subseq


            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < data.getDimen(); ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut.getT(0, i_prev)  + data.getCost(solut.getPos(i_prev), solut.getPos(j_next));
          cost_concat_2 = cost_concat_1 + solut.getT(j_next, k)  + data.getCost(solut.getPos(k), solut.getPos(i));
          cost_concat_3 = cost_concat_2 + solut.getT(i, j)       + data.getCost(solut.getPos(j), solut.getPos(k_next));

          cost_new = solut.getC(0, i_prev)                                                                  +   //       1st subseq
                  solut.getW(j_next, k)          * cost_concat_1 + solut.getC(j_next, k)             +   // concat 2nd subseq
                  solut.getW(i, j)               * cost_concat_2 + solut.getC(i, j)                  +   // concat 3rd subseq (reinserted seq)
                  solut.getW(k_next, data.getDimen()) * cost_concat_3 + solut.getC(k_next, data.getDimen());       // concat 4th subseq
          
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.getCost() - DBL_EPSILON) {
        solut.reinsert(I, J, POS+1);
        update_subseq_info_matrix(solut, data, I < POS+1 ? I : POS+1);

        solut.validate();

        if (solut.getCost() != cost_best) {
            std::cout << "deu merda\n";
        }

        return true;
    }

    return false;
}


void MLP::RVND(tSolution & solut, tData & data) {

    alignas(alignof(std::vector<int>)) std::vector<int> neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    alignas(INT_SZ) int index;
    alignas(INT_SZ) int neighbd;
    bool improve_flag;

    while (!neighbd_list.empty()) {

        index = data.getRndCrnt();
        neighbd = neighbd_list[index];

        improve_flag = false;

        switch(neighbd){
            case REINSERTION:
                improve_flag = search_reinsertion(solut, data, REINSERTION);
                break;				
            case OR_OPT_2:
                improve_flag = search_reinsertion(solut, data, OR_OPT_2);
                break;				
            case OR_OPT_3:
                improve_flag = search_reinsertion(solut, data, OR_OPT_3);
                break;				
            case SWAP:
                improve_flag = search_swap(solut, data);
                break;
            case TWO_OPT:
                improve_flag = search_two_opt(solut, data);
                break;				
        }
        if (improve_flag) {

            neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
        } else {
            
            neighbd_list.erase(neighbd_list.begin() + index);
        }
    }
}

std::vector<int> MLP::perturb(tSolution * solut, tData & data) {
    auto s = solut->getSolutVec();
    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        A_start = data.getRndCrnt();
        A_end = A_start + data.getRndCrnt();

        B_start = data.getRndCrnt();
        B_end = B_start + data.getRndCrnt();
    }

    auto reinsert = [](std::vector<int> &vec, int i, int j, int pos){
        std::vector<int> seq (vec.begin() + i, vec.begin() +j+1);
        if(pos < i){
            vec.erase(vec.begin() + i, vec.begin() + j+1);
            vec.insert(vec.begin() + pos, seq.begin(), seq.end());
        }else{
            vec.insert(vec.begin() + pos, seq.begin(), seq.end());
            vec.erase(vec.begin() + i, vec.begin() + j+1);
        }
    };
    
    if (A_start < B_start) {
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    } else {
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    return s;
}


void MLP::GILS_RVND(int Imax, int Iils, tData & data) {

    tSolution solut_partial(data);
    tSolution solut_crnt(data);
    tSolution solut_best(data);

    for(int i = 0; i < Imax; ++i){
        int aux = data.getRndCrnt();

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	

        solut_crnt.setSolutVec(
            construct(alpha, data)
        );

        update_subseq_info_matrix(solut_crnt, data);

        solut_partial.copy(solut_crnt);
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.getCost());	

        std::cout << "Solucao inicial ";
        solut_partial.print();

        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(solut_crnt, data);
            if(solut_crnt.getCost() < solut_partial.getCost() - DBL_EPSILON){
                solut_partial.copy(solut_crnt);
                iterILS = 0;

            }

            auto s = perturb(&solut_partial, data);
            solut_crnt.setSolutVec(s);
            update_subseq_info_matrix(solut_crnt, data);
            iterILS++;
        }

        if (solut_partial.getCost() < solut_best.getCost() - DBL_EPSILON) {
            solut_best.copy(solut_partial);
        }

        std::cout << "\tCurrent best cost: "<< solut_best.getCost() << std::endl;

        std::cout << "SOLUCAO: ";
        for(int i = 0; i < solut_best.getSolutVec().size(); i++){
            std::cout << solut_best.getPos(i) << " ";
        }
        std::cout << std::endl;

    }
    printf("COST: %.2lf\n", solut_best.getCost());
}
