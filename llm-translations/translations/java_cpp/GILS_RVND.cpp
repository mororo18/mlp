#include "GILS_RVND.h"
#include "Data.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

constexpr double GILS_RVND::R[];

GILS_RVND::GILS_RVND(bool verbose) {
    this->verbose = verbose;
    Data data = Data();
    data.loadData();

    int dimension = data.getDimension();
    std::vector<std::vector<double>> c(dimension, std::vector<double>(dimension, 0.0));
    for (int i = 0; i < dimension; i++) {
        for (int j = i; j < dimension; j++) {
            c[i][j] = data.getDistance(i, j);
            c[j][i] = data.getDistance(i, j);
        }
    }
    std::vector<int> & rnd = data.getRnd();
    info = std::make_unique<tInfo>(dimension, c, rnd);

    Iils = (dimension < 100 ? dimension : 100);
}

void GILS_RVND::subseq_load(tSolution & solut, tInfo & info) {
    int dimension = info.getDimen();
    for (int i = 0; i < dimension+1; i++) {
        int k = 1 - i - (i != 0 ? 0 : 1);

        solut.setSeq(i, i, info.T, 0.0);
        solut.setSeq(i, i, info.C, 0.0);
        solut.setSeq(i, i, info.W, (i != 0 ? 1.0 : 0.0));

        tArray & s = solut.getSolut();

        for (int j = i+1; j < dimension+1; j++) {
            int j_prev = j-1;

            solut.setSeq(i, j, info.T, info.getCost(s.get(j_prev), s.get(j)) + solut.getSeq(i, j_prev, info.T));
            solut.setSeq(i, j, info.C, solut.getSeq(i, j, info.T) + solut.getSeq(i, j_prev, info.C));
            solut.setSeq(i, j, info.W, j + k);
        }
    }

    solut.setCost(solut.getSeq(0, info.getDimen(), info.C));
}

tArray GILS_RVND::construction(double alpha, tInfo & info) {
    tArray s = tArray();
    s.reserve(info.getDimen()+1);
    s.add(0);

    tArray cList = tArray();
    cList.reserve(info.getDimen()-1);
    int dimension = info.getDimen();
    for (int i = 1; i < dimension; i++) {
        cList.add(i);
    }

    int r = 0;
    while (!cList.isEmpty()) {
        cList.sort(r, info);

        int index = info.rndCrnt();
        int c = cList.get(index);

        s.add(c);
        cList.remove(index);
        r = c;
    }

    s.add(0);
    return s;
}

bool GILS_RVND::search_swap(tSolution & solut, tInfo & info) {
    double cost_best = std::numeric_limits<double>::max();
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    double cost_concat_4;
    int I = -1;
    int J = -1;

    int dimension = info.getDimen();

    for (int i = 1; i < dimension-1; i++) {
        int i_prev = i - 1;
        int i_next = i + 1;
        int s_i_prev = solut.getSolut().get(i_prev);
        int s_i_next = solut.getSolut().get(i_next);
        int s_i = solut.getSolut().get(i);

        // immediate nodes case

        cost_concat_1 =                 solut.getSeq(0, i_prev, info.T) + info.getCost(s_i_prev, s_i_next);
        cost_concat_2 = cost_concat_1 + solut.getSeq(i, i_next, info.T) + info.getCost(s_i, solut.getSolut().get(i_next+1));

        cost_new = solut.getSeq(0, i_prev, info.C)                                                            // 1st subseq
                + solut.getSeq(i, i_next, info.W)           * (cost_concat_1) + info.getCost(s_i_next, s_i)              // concat 2nd subseq
                + solut.getSeq(i_next+1, dimension , info.W) * (cost_concat_2) + solut.getSeq(i_next+1, dimension, info.C);  // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - info.EPSILON;
            I = i;
            J = i_next;
        }

        for (int j = i_next+1; j < dimension; j++) {
            int j_next = j+1;
            int j_prev = j-1;
            int s_j = solut.getSolut().get(j);
            int s_j_prev = solut.getSolut().get(j_prev);
            int s_j_next = solut.getSolut().get(j_next);

            cost_concat_1 = solut.getSeq(0, i_prev, info.T) + info.getCost(s_i_prev, s_j);
            cost_concat_2 = cost_concat_1 + info.getCost(s_j, s_i_next);
            cost_concat_3 = cost_concat_2 + solut.getSeq(i_next, j_prev, info.T) + info.getCost(s_j_prev, s_i);
            cost_concat_4 = cost_concat_3 + info.getCost(s_i, s_j_next);

            cost_new = solut.getSeq(0, i_prev, info.C)                                                        /* first subseq */
                    + cost_concat_1                                                             /* concatenate second subseq (single node) */
                    + solut.getSeq(i_next, j_prev, info.W) * cost_concat_2 + solut.getSeq(i_next, j_prev, info.C)           /* concatenate third subseq */
                    + cost_concat_3                                                             /* concatenate fourth subseq (single node) */
                    + solut.getSeq(j_next, dimension, info.W) * cost_concat_4 + solut.getSeq(j_next, dimension, info.C);    /* concatenate fifth subseq */

            if (cost_new < cost_best) {
                cost_best = cost_new - info.EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.getCost()) {
        solut.getSolut().swap(I, J);
        subseq_load(solut, info);
        return true;
    }

    return false;
}

bool GILS_RVND::search_two_opt(tSolution & solut, tInfo & info) {
    int I = -1;
    int J = -1;
    double cost_best = std::numeric_limits<double>::max();
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    int dimension = info.getDimen();

    for (int i = 1; i < dimension-1; i++) {
        int i_prev = i -1;
        int s_i = solut.getSolut().get(i);
        int s_i_prev = solut.getSolut().get(i_prev);
        double rev_seq_cost = solut.getSeq(i, i+1, info.T);

        for (int j = i+2; j < dimension; j++) {
            int j_next = j+1;
            int s_j_next = solut.getSolut().get(j_next);
            int s_j_prev = solut.getSolut().get(j-1);
            int s_j = solut.getSolut().get(j);

            rev_seq_cost += info.getCost(s_j_prev, s_j) * (solut.getSeq(i, j, info.W)-1);

            cost_concat_1 =                 solut.getSeq(0, i_prev, info.T) + info.getCost(s_j, s_i_prev);
            cost_concat_2 = cost_concat_1 + solut.getSeq(i, j, info.T) + info.getCost(s_j_next, s_i);

            cost_new = solut.getSeq(0, i_prev, info.C)                                                        /*        1st subseq */
                    + solut.getSeq(i, j, info.W)              * cost_concat_1 + rev_seq_cost                  /* concat 2nd subseq (reversed seq) */
                    + solut.getSeq(j_next, dimension, info.W) * cost_concat_2 + solut.getSeq(j_next, dimension, info.C);    /* concat 3rd subseq */

            if (cost_new < cost_best) {
                cost_best = cost_new - info.EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.getCost() - info.EPSILON) {
        solut.getSolut().reverse(I, J);
        subseq_load(solut, info);
        return true;
    }

    return false;
}

bool GILS_RVND::search_reinsertion(tSolution & solut, tInfo & info, int opt) {
    double cost_best = std::numeric_limits<double>::max();
    double cost_new;
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    int I = -1;
    int J = -1;
    int POS = -1;
    int dimension = info.getDimen();

    for (int i = 1; i < dimension - opt + 1; i++) {
        int j = opt + i - 1;
        int i_prev = i-1;
        int j_next = j+1;
        int s_i = solut.getSolut().get(i);
        int s_j = solut.getSolut().get(j);
        int s_i_prev = solut.getSolut().get(i_prev);
        int s_j_next = solut.getSolut().get(j_next);

        // k -> reinsertion places
        for (int k = 0; k < i_prev; k++) {
            int k_next = k+1;
            int s_k = solut.getSolut().get(k);
            int s_k_next = solut.getSolut().get(k_next);

            cost_concat_1 = solut.getSeq(0, k, info.T) + info.getCost(s_k, s_i);
            cost_concat_2 = cost_concat_1 + solut.getSeq(i, j, info.T) + info.getCost(s_j, s_k_next);
            cost_concat_3 = cost_concat_2 + solut.getSeq(k_next, i_prev, info.T) + info.getCost(s_i_prev, s_j_next);

            cost_new = solut.getSeq(0, k, info.C)                                                             /*        1st subseq */
                    + solut.getSeq(i, j, info.W)              * cost_concat_1 + solut.getSeq(i, j, info.C)                  /* concat 2nd subseq (reinserted seq) */
                    + solut.getSeq(k_next, i_prev, info.W)    * cost_concat_2 + solut.getSeq(k_next, i_prev, info.C)        /* concat 3rd subseq */
                    + solut.getSeq(j_next, dimension, info.W) * cost_concat_3 + solut.getSeq(j_next, dimension, info.C);    /* concat 4th subseq */

            if (cost_new < cost_best) {
                cost_best = cost_new - info.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (int k = i+opt; k < dimension; k++) {
            int k_next = k+1;
            int s_k = solut.getSolut().get(k);
            int s_k_next = solut.getSolut().get(k_next);

            cost_concat_1 = solut.getSeq(0, i_prev, info.T) + info.getCost(s_i_prev, s_j_next);
            cost_concat_2 = cost_concat_1 + solut.getSeq(j_next, k, info.T) + info.getCost(s_k, s_i);
            cost_concat_3 = cost_concat_2 + solut.getSeq(i, j, info.T) + info.getCost(s_j, s_k_next);

            cost_new = solut.getSeq(0, i_prev, info.C)                                                        /*      1st subseq */
                    + solut.getSeq(j_next, k, info.W)         * cost_concat_1 + solut.getSeq(j_next, k, info.C)             /* concat 2nd subseq */
                    + solut.getSeq(i, j, info.W)              * cost_concat_2 + solut.getSeq(i, j, info.C)                  /* concat 3rd subseq (reinserted seq) */
                    + solut.getSeq(k_next, dimension, info.W) * cost_concat_3 + solut.getSeq(k_next, dimension, info.C);    /* concat 4th subseq */

            if (cost_new < cost_best) {
                cost_best = cost_new - info.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.getCost()) {
        solut.getSolut().reinsert(I, J, POS+1);
        subseq_load(solut, info);
        return true;
    }

    return false;
}

void GILS_RVND::RVND(tSolution & solut, tInfo & info) {
    tArray neighbd_list = tArray(std::vector<int>{
        info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT2, info.OR_OPT3
    });

    while (!neighbd_list.isEmpty()) {

        int i_rand = info.rndCrnt();
        int neighbd = neighbd_list.get(i_rand);

        bool improve = false;

        switch (neighbd) {
            case tInfo::REINSERTION:
                improve = search_reinsertion(solut, info, info.REINSERTION);
                break;
            case tInfo::OR_OPT2:
                improve = search_reinsertion(solut, info, info.OR_OPT2);
                break;
            case tInfo::OR_OPT3:
                improve = search_reinsertion(solut, info, info.OR_OPT3);
                break;
            case tInfo::SWAP:
                improve = search_swap(solut, info);
                break;
            case tInfo::TWO_OPT:
                improve = search_two_opt(solut, info);
                break;
        }

        if (improve) {
            neighbd_list.fill();
            neighbd_list.set(0, info.SWAP);
            neighbd_list.set(1, info.TWO_OPT);
            neighbd_list.set(2, info.REINSERTION);
            neighbd_list.set(3, info.OR_OPT2);
            neighbd_list.set(4, info.OR_OPT3);
        } else
            neighbd_list.remove(i_rand);
    }
}

tArray GILS_RVND::perturb(tArray & sl, tInfo & info) {
    tArray sl_cpy = sl.getArrayCopy();

    int A_start = 1, A_end = 1;
    int B_start = 1, B_end = 1;

    int size_max = (int) std::floor(sl.size_()/10);
    size_max = (size_max >= 2 ? size_max : 2);
    int size_min = 2;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = info.rndCrnt();
        A_end = A_start + info.rndCrnt();

        B_start = info.rndCrnt();
        B_end = B_start + info.rndCrnt();
    }

    if (A_start < B_start) {
        sl_cpy.reinsert(B_start, B_end-1, A_end);
        sl_cpy.reinsert(A_start, A_end-1, B_end);
    } else {
        sl_cpy.reinsert(A_start, A_end-1, B_end);
        sl_cpy.reinsert(B_start, B_end-1, A_end);
    }

    return sl_cpy;
}

void GILS_RVND::solve() {
    tInfo & info = *this->info;

    tSolution solut_best = tSolution(info.getDimen(), std::numeric_limits<double>::max());
    tSolution solut_crnt = tSolution(info.getDimen(), 0);
    tSolution solut_partial = tSolution(info.getDimen(), 0);

    for (int i = 0; i < Imax; i++) {
        int index = info.rndCrnt();

        double alpha = R[index];

        if (verbose) {
            std::cout << "[+] Local Search ";
            std::cout << i+1 << std::endl;
            std::cout << "\t[+] Constructing Inital Solution.." << std::endl;
        }

        solut_crnt.storeSolut(construction(alpha, info));
        subseq_load(solut_crnt, info);

        solut_partial.storeSolut(solut_crnt.getSolutCpy());
        solut_partial.setCost(solut_crnt.getCost());

        int iterILS = 0;
        if (verbose) {
            std::cout << "\t[+] Looking for the best Neighbor.." << std::endl;
        }
        while (iterILS < Iils) {
            RVND(solut_crnt, info);
            if (solut_crnt.getCost() < solut_partial.getCost()) {
                solut_partial.setCost(solut_crnt.getCost());
                solut_partial.storeSolut(solut_crnt.getSolutCpy());
                iterILS = 0;
            }

            solut_crnt.storeSolut(perturb(solut_partial.getSolut(), info));
            subseq_load(solut_crnt, info);
            iterILS++;
        }

        if (solut_partial.getCost() < solut_best.getCost()) {
            solut_best.storeSolut(solut_partial.getSolutCpy());
            solut_best.setCost(solut_partial.getCost());
        }

        if (verbose) {
            std::cout << "\tCurrent best solution cost: ";
            std::cout << std::setprecision(15) << solut_best.getCost() << std::endl;
        }
    }

    std::cout << "COST: ";
    std::cout << std::setprecision(15) << solut_best.getCost() << std::endl;
    std::cout << "SOLUTION: ";
    solut_best.getSolut().printPretty();
}
