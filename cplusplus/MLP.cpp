#include "MLP.hpp"

MLP::MLP(tInfo & info) {
    srand(clock());
    this->info = info;
}

void MLP::solve() {
    tSolution solut = Solution_init(info);

    int Imax = 10;
    int Iils;
    Iils = info.dimen < 100 ? info.dimen : 100;

    size_t start = clock();
    GILS_RVND(Imax, Iils, info);
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

std::vector<int> MLP::construct(const double alpha, tInfo & info){

    std::vector<int> s = {0};
    s.reserve(info.dimen+1);
    std::vector<int> cL(info.dimen-1);
    for(int i = 1; i < info.dimen; ++i){
        cL[i-1] = i;
    }

    int r = 0;
    while (!cL.empty()) {
        std::stable_sort(cL.begin(), cL.end(), 
            [r, &info] (const int i, const int j) {
                return info.cost[i][r] < info.cost[j][r];
            });

        /**/
        int range = std::ceil(cL.size() * alpha);
        int index = range > 0 ? rand() % range : 0;
        /**/

        //std::cout << info.rnd[info.rnd_index]<< std::endl;
        index = info.rnd[info.rnd_index++];
        int c = cL[index];
        s.push_back(c);
        //print_s(cL);
        r = c;
        cL.erase(cL.begin() + index);
    }

    s.push_back(0);

    return s;
}	

void MLP::subseq_load(tSolution & solut, const tInfo & info, int index = 0){
    alignas(INT_SZ) int i, j, j_prev, k;
    alignas(INT_SZ) int from = index;
    alignas(1) bool t;
    for (i = 0; i < info.dimen+1; i++) {
        k = 1 - i - (!i);
        t = i == from;

        solut.seq[i][i][info.T] = 0.0;
        solut.seq[i][i][info.C] = 0.0;
        solut.seq[i][i][info.W] = (double) !(i == 0);

        for (j = i+1; j < info.dimen+1; j++) {
            j_prev = j-1;
            
            solut.seq[i][j][info.T] = info.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][info.T];
            solut.seq[i][j][info.C] = solut.seq[i][j][info.T] + solut.seq[i][j_prev][info.C];
            solut.seq[i][j][info.W] = j + k;

        }
        from += t;
    }

    solut.cost = solut.seq[0][info.dimen][info.C];
}

bool MLP::search_swap(tSolution & solut, const tInfo & info) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    alignas(DBL_SZ) double cost_best = DBL_MAX;
    alignas(INT_SZ) int i, j, j_prev, j_next, i_prev, i_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < info.dimen-1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 solut.seq[0][i_prev][info.T] + info.cost[solut.s[i_prev]][solut.s[i_next]];
        cost_concat_2 = cost_concat_1 + solut.seq[i][i_next][info.T]  + info.cost[solut.s[i]][solut.s[i_next+1]];

        cost_new = solut.seq[0][i_prev][info.C]                                                    +           //       1st subseq
        solut.seq[i][i_next][info.W]               * (cost_concat_1) + info.cost[solut.s[i_next]][solut.s[i]]  +           // concat 2nd subseq
        solut.seq[i_next+1][info.dimen][info.W]   * (cost_concat_2) + solut.seq[i_next+1][info.dimen][info.C];   // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next+1; j < info.dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 solut.seq[0][i_prev][info.T]       + info.cost[solut.s[i_prev]][solut.s[j]];
            cost_concat_2 = cost_concat_1                           + info.cost[solut.s[j]][solut.s[i_next]];
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][j_prev][info.T]  + info.cost[solut.s[j_prev]][solut.s[i]];
            cost_concat_4 = cost_concat_3                           + info.cost[solut.s[i]][solut.s[j_next]];

            cost_new = solut.seq[0][i_prev][info.C]                                                 +      // 1st subseq
            cost_concat_1 +                                                             // concat 2nd subseq (single node)
            solut.seq[i_next][j_prev][info.W]      * cost_concat_2 + solut.seq[i_next][j_prev][info.C] +      // concat 3rd subseq
            cost_concat_3 +                                                             // concat 4th subseq (single node)
            solut.seq[j_next][info.dimen][info.W] * cost_concat_4 + solut.seq[j_next][info.dimen][info.C];   // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        swap(solut.s, I, J);
        subseq_load(solut, info, I);
        return true;
    }

    return false;
}

bool MLP::search_two_opt(tSolution & solut, const tInfo & info) {
    alignas(DBL_SZ) double cost_new, 
        cost_concat_1, cost_concat_2;
    alignas(DBL_SZ) double cost_best = DBL_MAX;// cost_l1, cost_l2;
    alignas(DBL_SZ) double rev_seq_cost;
    alignas(INT_SZ) int i, j, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;

    for (i = 1; i < info.dimen-1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = solut.seq[i][i+1][info.T];
        for (j = i + 2; j < info.dimen; ++j) {
            j_next = j + 1;


          rev_seq_cost += info.cost[solut.s[j-1]][solut.s[j]] * (solut.seq[i][j][info.W]-1.0);

          cost_concat_1 =                 solut.seq[0][i_prev][info.T]   + info.cost[solut.s[j]][solut.s[i_prev]];
          cost_concat_2 = cost_concat_1 + solut.seq[i][j][info.T]        + info.cost[solut.s[j_next]][solut.s[i]];

          cost_new = solut.seq[0][i_prev][info.C]                                                        +   //  1st subseq
              solut.seq[i][j][info.W]                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
              solut.seq[j_next][info.dimen][info.W] * cost_concat_2 + solut.seq[j_next][info.dimen][info.C];      // concat 3rd subseq

            
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        reverse(solut.s, I, J);
        subseq_load(solut, info);
        return true;
    }

    return false;
}

bool MLP::search_reinsertion(tSolution & solut, const tInfo & info, const int opt) {
    alignas(DBL_SZ) double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    alignas(DBL_SZ) double cost_best = DBL_MAX;//, cost_l1, cost_l2, cost_l3;
    alignas(INT_SZ) int i, j, k, k_next, i_prev, j_next;
    alignas(INT_SZ) int I;
    alignas(INT_SZ) int J;
    alignas(INT_SZ) int POS;

    for (i = 1, j = opt +i-1; i < info.dimen-opt+1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges 
        for(k = 0; k < i_prev; ++k){
            k_next = k+1;

          cost_concat_1 =                 solut.seq[0][k][info.T]            + info.cost[solut.s[k]][solut.s[i]];
          cost_concat_2 = cost_concat_1 + solut.seq[i][j][info.T]            + info.cost[solut.s[j]][solut.s[k_next]];
          cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev][info.T]  + info.cost[solut.s[i_prev]][solut.s[j_next]];

          cost_new = solut.seq[0][k][info.C]                                                                   +   //       1st subseq
              solut.seq[i][j][info.W]               * cost_concat_1 + solut.seq[i][j][info.C]                  +   //  concat 2nd subseq (reinserted seq)
              solut.seq[k_next][i_prev][info.W]     * cost_concat_2 + solut.seq[k_next][i_prev][info.C]        +   //  concat 3rd subseq
              solut.seq[j_next][info.dimen][info.W] * cost_concat_3 + solut.seq[j_next][info.dimen][info.C];       // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < info.dimen; ++k) {
            k_next = k + 1;

          cost_concat_1 =                 solut.seq[0][i_prev][info.T]  + info.cost[solut.s[i_prev]][solut.s[j_next]];
          cost_concat_2 = cost_concat_1 + solut.seq[j_next][k][info.T]  + info.cost[solut.s[k]][solut.s[i]];
          cost_concat_3 = cost_concat_2 + solut.seq[i][j][info.T]       + info.cost[solut.s[j]][solut.s[k_next]];

          cost_new = solut.seq[0][i_prev][info.C]                                                                  +   //       1st subseq
                  solut.seq[j_next][k][info.W]          * cost_concat_1 + solut.seq[j_next][k][info.C]             +   // concat 2nd subseq
                  solut.seq[i][j][info.W]               * cost_concat_2 + solut.seq[i][j][info.C]                  +   // concat 3rd subseq (reinserted seq)
                  solut.seq[k_next][info.dimen][info.W] * cost_concat_3 + solut.seq[k_next][info.dimen][info.C];       // concat 4th subseq
          
            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        reinsert(solut.s, I, J, POS+1);
        subseq_load(solut, info, I < POS+1 ? I : POS+1);
        return true;
    }

    return false;
}


void MLP::RVND(tSolution & solut, tInfo & info) {

    alignas(alignof(std::vector<int>)) std::vector<int> neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    alignas(INT_SZ) uint index;
    alignas(INT_SZ) int neighbd;
    //int k = 0;
    bool improve_flag;

    //cout << "RVND" << endl;
    while (!neighbd_list.empty()) {
        //k++;

        index = rand() % neighbd_list.size();
        //cout << info.rnd[info.rnd_index] << endl;
        index = info.rnd[info.rnd_index++];
        neighbd = neighbd_list[index];
        //std::cout <<"aq\n";

        improve_flag = false;

        switch(neighbd){
            case REINSERTION:
                //before();
                improve_flag = search_reinsertion(solut, info, REINSERTION);
                //after(REINSERTION);
                break;				
            case OR_OPT_2:
                //before();
                improve_flag = search_reinsertion(solut, info, OR_OPT_2);
                //after(OR_OPT2);
                break;				
            case OR_OPT_3:
                //before();
                improve_flag = search_reinsertion(solut, info, OR_OPT_3);
                //after(OR_OPT3);
                break;				
            case SWAP:
                //before();
                improve_flag = search_swap(solut, info);
                //after(SWAP);
                break;
            case TWO_OPT:
                //before();
                improve_flag = search_two_opt(solut, info);
                //after(TWO_OPT);
                break;				
        }
        //std::cout << (improve_flag ? "True" : "False") << std::endl;
        if (improve_flag) {
            neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
        } else {
            //std::cout << index << "  " << neighbd_list.size() << std::endl;
            //std::cout << solut.cost << std::endl;
            
            //std::cout << info.rnd_index << std::endl;
            neighbd_list.erase(neighbd_list.begin() + index);
        }

        //std::cout << "cost  " << solut.cost << std::endl ;


    }

    //exit(0);
    //std::cout << k << " RVND iteracoes" << std::endl;
}

std::vector<int> MLP::perturb(tSolution * solut, tInfo & info) {
    auto s = solut->s;
    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    int size_max = std::floor((info.dimen+1)/10);
    size_max = size_max >= 2 ? size_max : 2;
    int size_min = 2;
    //std::cout << "perturbing\n";
    //print_s(s);
    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {
        /**/
        int max = (info.dimen+1) -2 -size_max;
        A_start = rand() % max + 1;
        A_end = A_start + rand() % (size_max - size_min + 1) + size_min;

        B_start = rand() % max + 1;
        B_end = B_start + rand() % (size_max - size_min + 1) + size_min;
        /**/



        //std::cout << "paa\n";

        //cout << info.rnd[info.rnd_index] << endl;
        A_start = info.rnd[info.rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        A_end = A_start + info.rnd[info.rnd_index++];
        //std::cout << "A start  " << A_start << std::endl;
        //std::cout << "A end  " << A_end << std::endl;

        //cout << info.rnd[info.rnd_index] << endl;
        B_start = info.rnd[info.rnd_index++];
        //cout << info.rnd[info.rnd_index] << endl;
        B_end = B_start + info.rnd[info.rnd_index++];
        //std::cout << "B start  " << B_start << std::endl;
        //std::cout << "B end  " << B_end << std::endl;
    }
    
    //cout << "A_end  " << A_end << endl << "B_end  " << B_end << endl;

    if (A_start < B_start) {
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    } else {
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    //print_s(s);
    //subseq_load(solut, info);

    return s;
}


void MLP::GILS_RVND(int Imax, int Iils, tInfo & info) {

    tSolution solut_partial = Solution_init(info);
    tSolution solut_crnt = Solution_init(info);
    tSolution solut_best = Solution_init(info);

    for(int i = 0; i < Imax; ++i){
        /**/ int aux = (unsigned)rand() % TABLE_SZ;
        aux = info.rnd[info.rnd_index++];

        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	


        solut_crnt.s = construct(alpha, info);
        //print_s(solut_crnt.s);
        subseq_load(solut_crnt, info);

        //solut_partial = solut_crnt;
        Solution_cpy(solut_crnt, solut_partial, info);
        printf("\t[+] Looking for the best Neighbor..\n");
        printf("\t    Construction Cost: %.3lf\n", solut_partial.cost);	

        int iterILS = 0;
        //int k = 0;
        while (iterILS < Iils) {
            //k++;
            RVND(solut_crnt, info);
            if(solut_crnt.cost < solut_partial.cost - DBL_EPSILON){
                Solution_cpy(solut_crnt, solut_partial, info);
                //solut_partial = solut_crnt;
                iterILS = 0;
            }

            solut_crnt.s = perturb(&solut_partial, info);
            subseq_load(solut_crnt, info);
            //exit(0);
            //std::cout << "ITER  " << iterILS << std::endl;
            iterILS++;
        }

        //subseq_load(solut_partial, info);

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(solut_partial, solut_best, info);
            //solut_best = solut_partial;
        }

        //after(7);

        //std::cout << "\tCurrent search cost: "<< cost_sl << std::endl;
        std::cout << "\tCurrent best cost: "<< solut_best.cost << std::endl;
        //std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        //std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;
        //std::cout << k << "  Iteracoes " << std::endl;

        std::cout << "SOLUCAO: ";
        for(int i = 0; i < solut_best.s.size(); i++){
            std::cout << solut_best.s[i] << " ";
        }
        std::cout << std::endl;

    }
    //std::cout << "Dimension: " << dimension << std::endl;
    printf("COST: %.2lf\n", solut_best.cost);
}

