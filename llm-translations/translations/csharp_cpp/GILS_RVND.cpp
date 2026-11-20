#include "GILS_RVND.hpp"
#include "tInfo.hpp"
#include "dbl.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace MLP {

    GILS_RVND::GILS_RVND() {
        data.loadData();

        std::cout << "dimesion: " << data.GetDimen() << "." << std::endl;

        int dimension = data.GetDimen();

        Iils = (dimension < 100 ? dimension : 100);

        const std::vector<int> & rnd = data.GetRnd();
        (void)rnd;
    }

    void GILS_RVND::update_subseq_info_matrix(tSolution & solut) {
        std::vector<int> & solutSolut = solut.GetSolut();
        for (int i = 0; i < data.GetDimen() + 1; i++) {
            int k = 1 - i - (i != 0 ? 0 : 1);

            solut.SetSeq(i, i, tInfoIndex::T, 0.0);
            solut.SetSeq(i, i, tInfoIndex::C, 0.0);
            solut.SetSeq(i, i, tInfoIndex::W, (i != 0 ? 1.0 : 0.0));

            for (int j = i + 1; j < data.GetDimen() + 1; j++) {
                int j_prev = j - 1;

                solut.SetSeq(i, j, tInfoIndex::T, data.GetCost(solutSolut[j_prev], solutSolut[j]) + solut.GetSeq(i, j_prev, tInfoIndex::T));
                solut.SetSeq(i, j, tInfoIndex::C, solut.GetSeq(i, j, tInfoIndex::T) + solut.GetSeq(i, j_prev, tInfoIndex::C));
                solut.SetSeq(i, j, tInfoIndex::W, j + k);
            }
        }

        solut.SetCost(solut.GetSeq(0, data.GetDimen(), tInfoIndex::C));
    }

    void GILS_RVND::sort(std::vector<int> & arr, int r) {
        Quicksort(arr, 0, (int)arr.size() - 1, r);
    }

    void GILS_RVND::Quicksort(std::vector<int> & arr, int left, int right, int r) {
        if (left < right) {
            int pivot = Partition(arr, left, right, r);
            Quicksort(arr, left, pivot - 1, r);
            Quicksort(arr, pivot + 1, right, r);
        }
    }

    int GILS_RVND::Partition(std::vector<int> & arr, int left, int right, int r) {
        int pivotIndex = right;
        double pivotValue = data.GetCost(r, arr[pivotIndex]);
        int i = left - 1;

        for (int j = left; j < right; j++) {
            if (data.GetCost(r, arr[j]) < pivotValue) {
                i++;
                std::swap(arr[i], arr[j]); // Troca os elementos
            }
        }
        std::swap(arr[i + 1], arr[right]); // Troca o pivo
        return i + 1;
    }

    std::vector<int> GILS_RVND::construction(double alpha) {
        (void)alpha;
        std::vector<int> s = {0};

        std::vector<int> cList;
        for (int i = 1; i < data.GetDimen(); i++) {
            cList.push_back(i);
        }

        int r = 0;
        while (cList.size() > 0) {

            sort(cList, r);

            int r_value = data.GetRndCrnt();

            int cN = cList[r_value];
            s.push_back(cN);
            cList.erase(std::find(cList.begin(), cList.end(), cN));
            r = cN;
        }
        s.push_back(0);

        return s;
    }

    void GILS_RVND::swap(std::vector<int> & s, int i, int j) {
        int tmp = s[i];
        s[i] = s[j];
        s[j] = tmp;
    }

    void GILS_RVND::reverse(std::vector<int> & s, int i, int j) {
        std::reverse(s.begin() + i, s.begin() + j + 1);
    }

    void GILS_RVND::reinsert(std::vector<int> & s, int i, int j, int pos) {
        int sz = j - i + 1;
        std::vector<int> range(s.begin() + i, s.begin() + i + sz);
        if (i < pos) {
            s.insert(s.begin() + pos, range.begin(), range.end());
            s.erase(s.begin() + i, s.begin() + i + sz);
        } else {
            s.insert(s.begin() + pos, range.begin(), range.end());
            s.erase(s.begin() + i + sz, s.begin() + i + sz + sz);
        }
    }

    bool GILS_RVND::search_swap(tSolution & solut) {
        double cost_best = std::numeric_limits<double>::max();
        double cost_new;
        double cost_concat_1, cost_concat_2,
               cost_concat_3, cost_concat_4;
        int I = -1;
        int J = -1;

        std::vector<int> & solutSolut = solut.GetSolut();
        int dimen = data.GetDimen();

        for (int i = 1; i < dimen - 1; i++) {
            int i_prev = i - 1;
            int i_next = i + 1;

            cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex::T) + data.GetCost(solut.GetPos(i_prev), solut.GetPos(i_next));
            cost_concat_2 = cost_concat_1 + solut.GetSeq(i, i_next, tInfoIndex::T) + data.GetCost(solutSolut[i], solutSolut[i_next + 1]);

            cost_new = solut.GetSeq(0, i_prev, tInfoIndex::C)
                    + solut.GetSeq(i, i_next, tInfoIndex::W)         * (cost_concat_1) + data.GetCost(solutSolut[i_next], solutSolut[i])
                    + solut.GetSeq(i_next + 1, dimen, tInfoIndex::W) * (cost_concat_2) + solut.GetSeq(i_next + 1, dimen, tInfoIndex::C);

            if (cost_new < cost_best) {
                cost_best = cost_new - EPSILON;
                I = i;
                J = i_next;
            }

            for (int j = i_next + 1; j < dimen; j++) {
                int j_next = j + 1;
                int j_prev = j - 1;

                cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex::T) + data.GetCost(solutSolut[i_prev], solutSolut[j]);
                cost_concat_2 = cost_concat_1 + data.GetCost(solutSolut[j], solutSolut[i_next]);
                cost_concat_3 = cost_concat_2 + solut.GetSeq(i_next, j_prev, tInfoIndex::T) + data.GetCost(solutSolut[j_prev], solutSolut[i]);
                cost_concat_4 = cost_concat_3 + data.GetCost(solutSolut[i], solutSolut[j_next]);

                cost_new = solut.GetSeq(0, i_prev, tInfoIndex::C)
                        + cost_concat_1
                        + solut.GetSeq(i_next, j_prev, tInfoIndex::W) * cost_concat_2 + solut.GetSeq(i_next, j_prev, tInfoIndex::C)
                        + cost_concat_3
                        + solut.GetSeq(j_next, dimen, tInfoIndex::W) * cost_concat_4 + solut.GetSeq(j_next, dimen, tInfoIndex::C);

                if (cost_new < cost_best) {
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                }
            }
        }

        if (cost_best < solut.GetCost() - EPSILON) {
            swap(solutSolut, I, J);
            update_subseq_info_matrix(solut);
            return true;
        }

        return false;
    }

    bool GILS_RVND::search_two_opt(tSolution & solut) {
        double cost_best = std::numeric_limits<double>::max();
        double cost_new;
        double cost_concat_1,
               cost_concat_2;
        int I = -1;
        int J = -1;
        std::vector<int> & solutSolut = solut.GetSolut();

        int dimen = data.GetDimen();

        for (int i = 1; i < dimen - 1; i++) {
            int i_prev = i - 1;

            double rev_seq_cost = solut.GetSeq(i, i + 1, tInfoIndex::T);

            for (int j = i + 2; j < dimen; j++) {
                int j_next = j + 1;

                rev_seq_cost += data.GetCost(solutSolut[j - 1], solutSolut[j]) * (solut.GetSeq(i, j, tInfoIndex::W) - 1.0);

                cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex::T) + data.GetCost(solutSolut[j], solutSolut[i_prev]);
                cost_concat_2 = cost_concat_1 + solut.GetSeq(i, j, tInfoIndex::T) + data.GetCost(solutSolut[j_next], solutSolut[i]);

                cost_new = solut.GetSeq(0, i_prev, tInfoIndex::C)
                        + solut.GetSeq(i, j, tInfoIndex::W)         * cost_concat_1 + rev_seq_cost
                        + solut.GetSeq(j_next, dimen, tInfoIndex::W) * cost_concat_2 + solut.GetSeq(j_next, dimen, tInfoIndex::C);

                if (cost_new < cost_best) {
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                }
            }
        }

        if (cost_best < solut.GetCost() - EPSILON) {
            reverse(solutSolut, I, J);
            update_subseq_info_matrix(solut);
            return true;
        }

        return false;
    }

    bool GILS_RVND::search_reinsertion(tSolution & solut, int opt) {
        double cost_best = std::numeric_limits<double>::max();
        double cost_new;
        double cost_concat_1, cost_concat_2,
               cost_concat_3;
        int I = -1;
        int J = -1;
        int POS = -1;
        std::vector<int> & solutSolut = solut.GetSolut();

        int dimen = data.GetDimen();

        for (int i = 1; i < dimen - opt + 1; i++) {
            int j = opt + i - 1;
            int i_prev = i - 1;
            int j_next = j + 1;

            for (int k = 0; k < i_prev; k++) {
                int k_next = k + 1;

                cost_concat_1 = solut.GetSeq(0, k, tInfoIndex::T) + data.GetCost(solutSolut[k], solutSolut[i]);
                cost_concat_2 = cost_concat_1 + solut.GetSeq(i, j, tInfoIndex::T) + data.GetCost(solutSolut[j], solutSolut[k_next]);
                cost_concat_3 = cost_concat_2 + solut.GetSeq(k_next, i_prev, tInfoIndex::T) + data.GetCost(solutSolut[i_prev], solutSolut[j_next]);

                cost_new = solut.GetSeq(0, k, tInfoIndex::C)
                        + solut.GetSeq(i, j, tInfoIndex::W)           * cost_concat_1 + solut.GetSeq(i, j, tInfoIndex::C)
                        + solut.GetSeq(k_next, i_prev, tInfoIndex::W) * cost_concat_2 + solut.GetSeq(k_next, i_prev, tInfoIndex::C)
                        + solut.GetSeq(j_next, dimen, tInfoIndex::W)  * cost_concat_3 + solut.GetSeq(j_next, dimen, tInfoIndex::C);

                if (cost_new < cost_best) {
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                    POS = k;
                }
            }

            for (int k = i + opt; k < dimen; k++) {
                int k_next = k + 1;

                cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex::T) + data.GetCost(solutSolut[i_prev], solutSolut[j_next]);
                cost_concat_2 = cost_concat_1 + solut.GetSeq(j_next, k, tInfoIndex::T) + data.GetCost(solutSolut[k], solutSolut[i]);
                cost_concat_3 = cost_concat_2 + solut.GetSeq(i, j, tInfoIndex::T) + data.GetCost(solutSolut[j], solutSolut[k_next]);

                cost_new = solut.GetSeq(0, i_prev, tInfoIndex::C)
                        + solut.GetSeq(j_next, k, tInfoIndex::W)     * cost_concat_1 + solut.GetSeq(j_next, k, tInfoIndex::C)
                        + solut.GetSeq(i, j, tInfoIndex::W)          * cost_concat_2 + solut.GetSeq(i, j, tInfoIndex::C)
                        + solut.GetSeq(k_next, dimen, tInfoIndex::W) * cost_concat_3 + solut.GetSeq(k_next, dimen, tInfoIndex::C);

                if (cost_new < cost_best) {
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                    POS = k;
                }
            }
        }

        if (cost_best < solut.GetCost() - EPSILON) {
            reinsert(solutSolut, I, J, POS + 1);
            update_subseq_info_matrix(solut);
            return true;
        }

        return false;
    }

    void GILS_RVND::RVND(tSolution & solut) {
        std::vector<int> neighbd_list = {Neighborhood::SWAP, Neighborhood::TWO_OPT, Neighborhood::REINSERTION, Neighborhood::OR_OPT2, Neighborhood::OR_OPT3};

        /*
        var t = new List<int>();
        for(int i = 0; i < dimension; i++)
            t.Add(i);
        t.Add(0);
        */

        while (neighbd_list.size() != 0) {
            int i_rand = data.GetRndCrnt();

            int neighbd = neighbd_list[i_rand];

            bool improve = false;

            switch (neighbd) {
                case Neighborhood::REINSERTION:
                    improve = search_reinsertion(solut, Neighborhood::REINSERTION);
                    break;
                case Neighborhood::OR_OPT2:
                    improve = search_reinsertion(solut, Neighborhood::OR_OPT2);
                    break;
                case Neighborhood::OR_OPT3:
                    improve = search_reinsertion(solut, Neighborhood::OR_OPT3);
                    break;
                case Neighborhood::SWAP:
                    improve = search_swap(solut);
                    break;
                case Neighborhood::TWO_OPT:
                    improve = search_two_opt(solut);
                    break;
            }

            if (improve) {
                neighbd_list = {Neighborhood::SWAP, Neighborhood::TWO_OPT, Neighborhood::REINSERTION, Neighborhood::OR_OPT2, Neighborhood::OR_OPT3};
            } else {
                neighbd_list.erase(neighbd_list.begin() + i_rand);
            }
        }
    }

    std::vector<int> GILS_RVND::perturb(const std::vector<int> & sl) {
        std::vector<int> s(sl);

        int A_start = 1, A_end = 1;
        int B_start = 1, B_end = 1;

        int size_max = (int)std::floor((double)sl.size() / 10);
        size_max = (size_max >= 2 ? size_max : 2);
        int size_min = 2;
        (void)size_max;
        (void)size_min;

        while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

            A_start = data.GetRndCrnt();
            A_end = A_start + data.GetRndCrnt();

            B_start = data.GetRndCrnt();
            B_end = B_start + data.GetRndCrnt();
        }

        if (A_start < B_start) {
            reinsert(s, B_start, B_end - 1, A_end);
            reinsert(s, A_start, A_end - 1, B_end);
        } else {
            reinsert(s, A_start, A_end - 1, B_end);
            reinsert(s, B_start, B_end - 1, A_end);
        }

        return s;
    }

    void GILS_RVND::solve(bool verbose) {
        tSolution solut_best(data.GetDimen(), std::numeric_limits<double>::max());
        tSolution solut_crnt(data.GetDimen(), 0.0);
        tSolution solut_partial(data.GetDimen(), 0.0);

        for (int i = 0; i < Imax; i++) {
            int index = data.GetRndCrnt();
            double alpha = R[index];

            if (verbose) {
                std::cout << "[+] Local Search " << (i + 1) << std::endl;
                std::cout << "\t[+] Constructing Inital Solution.." << std::endl;
            }

            solut_crnt.StoreSolut(construction(alpha));

            update_subseq_info_matrix(solut_crnt);

            solut_partial.StoreSolut(solut_crnt.GetSolutCpy());
            solut_partial.SetCost(solut_crnt.GetCost());

            if (verbose) {
                std::cout << "Construction: " << DblToStr(solut_crnt.GetCost()) << std::endl;
                std::cout << "\t[+] Looking for the best Neighbor.." << std::endl;
            }
            int iterILS = 0;
            while (iterILS < Iils) {
                RVND(solut_crnt);
                if (solut_crnt.GetCost() < solut_partial.GetCost()) {
                    solut_partial.StoreSolut(solut_crnt.GetSolutCpy());
                    solut_partial.SetCost(solut_crnt.GetCost());

                    iterILS = 0;
                }

                solut_crnt.StoreSolut(perturb(solut_partial.GetSolut()));
                update_subseq_info_matrix(solut_crnt);
                iterILS++;
            }

            if (solut_partial.GetCost() < solut_best.GetCost()) {
                solut_best.StoreSolut(solut_partial.GetSolutCpy());
                solut_best.SetCost(solut_partial.GetCost());
            }

            if (verbose) {
                std::cout << "\tCurrent best solution cost: " << DblToStr(solut_best.GetCost()) << std::endl;
            }
        }

        std::cout << "SOLUTION: (";
        const std::vector<int> & best = solut_best.GetSolut();
        for (size_t k = 0; k < best.size(); k++) {
            if (k != 0)
                std::cout << ", ";
            std::cout << best[k];
        }
        std::cout << ")." << std::endl;
        std::cout << "COST: " << DblToStr(solut_best.GetCost()) << std::endl;
    }

}
