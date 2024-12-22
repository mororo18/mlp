using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;

namespace MLP {
    class GILS_RVND {

        public const double EPSILON = 1e-15;

        Data data;

        private  int                    Iils;
        private const int               Imax = 10;
        private double []  R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
        private const int               R_size = 26;

        public GILS_RVND(){
            data = new Data();
            data.loadData();

            Console.WriteLine(string.Format("dimesion: {0}.", data.GetDimen()));

            int dimension = data.GetDimen();

            Iils = (dimension < 100 ? dimension : 100);


            int [] rnd = data.GetRnd();
            
    
        }

        private void subseq_load(tSolution solut){
            var solutSolut = solut.GetSolut();  
            for(int i = 0; i < data.GetDimen()+1; i++){
                int k = 1 - i - (i != 0 ? 0 : 1);

                solut.SetSeq(i, i, tInfoIndex.T, 0.0);
                solut.SetSeq(i, i, tInfoIndex.C, 0.0);
                solut.SetSeq(i, i, tInfoIndex.W, (i != 0 ? 1.0 : 0.0));

                for(int j = i+1; j < data.GetDimen()+1; j++){
                    int j_prev = j-1;

                    solut.SetSeq(i, j, tInfoIndex.T, data.GetCost(solutSolut[j_prev], solutSolut[j]) + solut.GetSeq(i, j_prev, tInfoIndex.T));
                    solut.SetSeq(i, j, tInfoIndex.C, solut.GetSeq(i, j, tInfoIndex.T) + solut.GetSeq(i, j_prev, tInfoIndex.C));
                    solut.SetSeq(i, j, tInfoIndex.W, j + k);

                }
            }

            solut.SetCost(solut.GetSeq(0, data.GetDimen(), tInfoIndex.C));
        }


        private void sort(List<int> arr, int r)
        {
            Quicksort(arr, 0, arr.Count - 1, r);
        }

        private void Quicksort(List<int> arr, int left, int right, int r)
        {
            if (left < right)
            {
                int pivot = Partition(arr, left, right, r);
                Quicksort(arr, left, pivot - 1, r);
                Quicksort(arr, pivot + 1, right, r);
            }
        }

        private int Partition(List<int> arr, int left, int right, int r)
        {
            int pivotIndex = right;
            double pivotValue = data.GetCost(r, arr[pivotIndex]);
            int i = left - 1;

            for (int j = left; j < right; j++)
            {
                if (data.GetCost(r, arr[j]) < pivotValue)
                {
                    i++;
                    (arr[i], arr[j]) = (arr[j], arr[i]); // Troca os elementos usando tuple swap
                }
            }
            (arr[i + 1], arr[right]) = (arr[right], arr[i + 1]); // Troca o pivô
            return i + 1;
        }

        

        private List<int> construction(double alpha){
            var s = new List<int> {0};

            var cList = new List<int>();
            for(int i = 1; i < data.GetDimen(); i++){
                cList.Add(i);
            }

            int r = 0;
            while(cList.Count > 0){

                sort(cList, r);

                int r_value = data.GetRndCrnt();

                int cN = cList[r_value];
                s.Add(cN);
                cList.Remove(cN);
                r = cN;
            }
            s.Add(0);

            return s;
        }

        private void swap(List<int> s, int i, int j){
            //Console.WriteLine(string.Format("Swap: ({0}, {1}).", i, j));
            int tmp = s[i];
            s[i] = s[j];
            s[j] = tmp;
            //Console.WriteLine(string.Format("Swap: ({0}).", string.Join(", ", s)));
        }

        private void reverse(List<int> s, int i, int j){
            s.Reverse(i, j-i+1);
            //Console.WriteLine(string.Format("Reverse: ({0}).", string.Join(", ", s)));
        }

        private void reinsert(List<int> s, int i, int j, int pos){
            int sz = j-i+1;
            if (i < pos) {
                s.InsertRange(pos, s.GetRange(i, sz));
                s.RemoveRange(i, sz);
            } else {
                s.InsertRange(pos, s.GetRange(i, sz));
                s.RemoveRange(i+sz, sz);
            }
            //Console.WriteLine(string.Format("Reinsert: ({0}).", string.Join(", ", s)));
        }

        private bool search_swap(tSolution solut){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1, cost_concat_2,
                   cost_concat_3, cost_concat_4;
            int I = -1;
            int J = -1;

            var solutSolut = solut.GetSolut();
            var dimen = data.GetDimen();

            for(int i = 1; i < dimen-1; i++){
                int i_prev = i - 1;
                int i_next = i + 1;

                cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex.T) + data.GetCost(solut.GetPos(i_prev),  solut.GetPos(i_next));
                cost_concat_2 = cost_concat_1 + solut.GetSeq(i, i_next, tInfoIndex.T) + data.GetCost(solutSolut[i], solutSolut[i_next+1]);

                cost_new = solut.GetSeq(0, i_prev, tInfoIndex.C)
                        + solut.GetSeq(i, i_next, tInfoIndex.W)             * (cost_concat_1) + data.GetCost(solutSolut[i_next], solutSolut[i])
                        + solut.GetSeq(i_next+1, dimen, tInfoIndex.W)   * (cost_concat_2) + solut.GetSeq(i_next+1, dimen, tInfoIndex.C);

                if(cost_new < cost_best){
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = i_next;
                }

                for(int j = i_next+1; j < dimen; j++){
                    int j_next = j+1;
                    int j_prev = j-1;

                    cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex.T) + data.GetCost(solutSolut[i_prev], solutSolut[j]);
                    cost_concat_2 = cost_concat_1 + data.GetCost(solutSolut[j], solutSolut[i_next]);
                    cost_concat_3 = cost_concat_2 + solut.GetSeq(i_next, j_prev, tInfoIndex.T) + data.GetCost(solutSolut[j_prev], solutSolut[i]);
                    cost_concat_4 = cost_concat_3 + data.GetCost(solutSolut[i], solutSolut[j_next]);

                    cost_new = solut.GetSeq(0, i_prev, tInfoIndex.C)
                            + cost_concat_1
                            + solut.GetSeq(i_next, j_prev, tInfoIndex.W) * cost_concat_2 + solut.GetSeq(i_next, j_prev, tInfoIndex.C)
                            + cost_concat_3
                            + solut.GetSeq(j_next, dimen, tInfoIndex.W) * cost_concat_4 + solut.GetSeq(j_next, dimen, tInfoIndex.C);
                    
                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if(cost_best < solut.GetCost() - EPSILON){
              //if (I == -1 && J == -1) {
              //    Console.WriteLine(string.Format("Swap: ({0}, {1}, {2}).", cost_best, solut.GetCost(), -cost_best+solut.GetCost()));
              //}
                swap(solutSolut, I, J);
                subseq_load(solut);
                return true;
            }

            return false;
        }

        private bool search_two_opt(tSolution solut){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1,
                   cost_concat_2;
            int I = -1;
            int J = -1;
            var solutSolut = solut.GetSolut();

            var dimen = data.GetDimen();

            for (int i = 1; i < dimen-1; i++){
                int i_prev = i-1;

                double rev_seq_cost = solut.GetSeq(i, i+1, tInfoIndex.T);

                for(int j = i+2; j < dimen; j++){
                    int j_next = j+1;
                    //int j_prev = j-1;

                    rev_seq_cost += data.GetCost(solutSolut[j-1], solutSolut[j]) * (solut.GetSeq(i,  j, tInfoIndex.W)-1.0);

                    cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex.T) + data.GetCost(solutSolut[j], solutSolut[i_prev]);
                    cost_concat_2 = cost_concat_1 + solut.GetSeq(i, j,  tInfoIndex.T) + data.GetCost(solutSolut[j_next], solutSolut[i]);

                    cost_new = solut.GetSeq(0, i_prev,  tInfoIndex.C)
                            + solut.GetSeq(i, j, tInfoIndex.W)              * cost_concat_1 + rev_seq_cost
                            + solut.GetSeq(j_next, dimen, tInfoIndex.W) * cost_concat_2 + solut.GetSeq(j_next, dimen, tInfoIndex.C);

                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i;
                        J = j;
                    }

                    //Environment.Exit(0);
                }
            }

            if(cost_best < solut.GetCost() - EPSILON){
                reverse(solutSolut, I, J);
                subseq_load(solut);
                return true;
            }

            return false;
        }

        private bool search_reinsertion(tSolution solut, int opt){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1, cost_concat_2,
                   cost_concat_3;
            int I = -1;
            int J = -1;
            int POS = -1;
            var solutSolut = solut.GetSolut();

            var dimen = data.GetDimen();

            for (int i = 1; i < dimen - opt + 1; i++){
                int j = opt + i - 1;
                int i_prev = i - 1;
                int j_next = j + 1;

                for(int k = 0; k < i_prev; k++){
                    int k_next = k+1;

                    cost_concat_1 = solut.GetSeq(0, k, tInfoIndex.T) + data.GetCost(solutSolut[k], solutSolut[i]);
                    cost_concat_2 = cost_concat_1 + solut.GetSeq(i, j, tInfoIndex.T) + data.GetCost(solutSolut[j], solutSolut[k_next]);
                    cost_concat_3 = cost_concat_2 + solut.GetSeq(k_next, i_prev, tInfoIndex.T) + data.GetCost(solutSolut[i_prev], solutSolut[j_next]);

                    cost_new = solut.GetSeq(0, k, tInfoIndex.C)
                            + solut.GetSeq(i, j, tInfoIndex.W)              * cost_concat_1 + solut.GetSeq(i, j, tInfoIndex.C)
                            + solut.GetSeq(k_next, i_prev, tInfoIndex.W)    * cost_concat_2 + solut.GetSeq(k_next, i_prev, tInfoIndex.C)
                            + solut.GetSeq(j_next, dimen, tInfoIndex.W) * cost_concat_3 + solut.GetSeq(j_next, dimen, tInfoIndex.C);

                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }

                for(int k = i+opt; k < dimen; k++){
                    int k_next = k+1;

                    cost_concat_1 = solut.GetSeq(0, i_prev, tInfoIndex.T) + data.GetCost(solutSolut[i_prev], solutSolut[j_next]);
                    cost_concat_2 = cost_concat_1 + solut.GetSeq(j_next, k, tInfoIndex.T) + data.GetCost(solutSolut[k], solutSolut[i]);
                    cost_concat_3 = cost_concat_2 + solut.GetSeq(i, j, tInfoIndex.T) + data.GetCost(solutSolut[j], solutSolut[k_next]);

                    cost_new = solut.GetSeq(0, i_prev, tInfoIndex.C)
                            + solut.GetSeq(j_next, k, tInfoIndex.W)         * cost_concat_1 + solut.GetSeq(j_next, k, tInfoIndex.C)
                            + solut.GetSeq(i, j, tInfoIndex.W)              * cost_concat_2 + solut.GetSeq(i, j, tInfoIndex.C)
                            + solut.GetSeq(k_next, dimen, tInfoIndex.W) * cost_concat_3 + solut.GetSeq(k_next, dimen, tInfoIndex.C);

                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i; 
                        J = j;
                        POS = k;
                    }
                }
            }

            if(cost_best < solut.GetCost() - EPSILON){
                //Console.WriteLine("Reinsertion");
                //Console.WriteLine(cost_best);
                reinsert(solutSolut, I, J, POS+1);
                subseq_load(solut);
                //Console.WriteLine(seq[0, dimension, C]);
                return true;
            }

            return false;

        }

        private void RVND(tSolution solut){
            List<int> neighbd_list = new List<int> {Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION, Neighborhood.OR_OPT2, Neighborhood.OR_OPT3};

            /*
            var t = new List<int>();
            for(int i = 0; i < dimension; i++)
                t.Add(i);
            t.Add(0);
            */

            while(neighbd_list.Count != 0){
                int i_rand = data.GetRndCrnt();

                int neighbd = neighbd_list[i_rand];

                bool improve = false;

                switch(neighbd){
                    case Neighborhood.REINSERTION:
                        improve = search_reinsertion(solut, Neighborhood.REINSERTION);
                        break;
                    case Neighborhood.OR_OPT2:
                        improve = search_reinsertion(solut, Neighborhood.OR_OPT2);
                        break;
                    case Neighborhood.OR_OPT3:
                        improve = search_reinsertion(solut, Neighborhood.OR_OPT3);
                        break;
                    case Neighborhood.SWAP:
                        improve = search_swap(solut);
                        break;
                    case Neighborhood.TWO_OPT:
                        improve = search_two_opt(solut);
                        break;
                }

                if(improve){
                    neighbd_list = new List<int> {Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION, Neighborhood.OR_OPT2, Neighborhood.OR_OPT3};
                }else{
                    neighbd_list.RemoveAt(i_rand);
                }
            }
        }

        private List<int> perturb(List<int> sl){
            var s = new List<int>(sl);

            int A_start = 1, A_end = 1;
            int B_start = 1, B_end = 1;

            int size_max = (int)Math.Floor((double)sl.Count/10);
            size_max = (size_max >= 2 ? size_max : 2);
            int size_min = 2;

            while((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)){

                A_start = data.GetRndCrnt();
                A_end = A_start + data.GetRndCrnt();
                
                B_start = data.GetRndCrnt();
                B_end = B_start + data.GetRndCrnt();
            }

            if(A_start < B_start) {
                reinsert(s, B_start, B_end-1, A_end);
                reinsert(s, A_start, A_end-1, B_end);
            } else {
                reinsert(s, A_start, A_end-1, B_end);
                reinsert(s, B_start, B_end-1, A_end);
            }

            return s;
        }


        public void solve(){
            var solut_best = new tSolution(data.GetDimen(), Double.MaxValue);
            var solut_crnt = new tSolution(data.GetDimen(), 0.0);
            var solut_partial = new tSolution(data.GetDimen(), 0.0);

            for(int i = 0; i < Imax; i++){
                int index = data.GetRndCrnt();
                double alpha = R[index];

                Console.WriteLine("[+] Local Search " + (i+1));
                Console.WriteLine("\t[+] Constructing Inital Solution..");

                solut_crnt.StoreSolut(construction(alpha));

                subseq_load(solut_crnt);

                solut_partial.StoreSolut(solut_crnt.GetSolutCpy());
                solut_partial.SetCost(solut_crnt.GetCost());

                Console.WriteLine("Construction: " + solut_crnt.GetCost());

                Console.WriteLine("\t[+] Looking for the best Neighbor..");
                int iterILS = 0;
                while(iterILS < Iils){
                    RVND(solut_crnt);
                    if(solut_crnt.GetCost() < solut_partial.GetCost()){
                        solut_partial.StoreSolut(solut_crnt.GetSolutCpy());
                        solut_partial.SetCost(solut_crnt.GetCost());

                        iterILS = 0;
                    }

                    solut_crnt.StoreSolut(perturb(solut_partial.GetSolut()));
                    subseq_load(solut_crnt);
                    iterILS++;
                }

                if(solut_partial.GetCost() < solut_best.GetCost()){
                    solut_best.StoreSolut(solut_partial.GetSolutCpy());
                    solut_best.SetCost(solut_partial.GetCost());
                }

                Console.WriteLine("\tCurrent best solution cost: "+solut_best.GetCost());
            }
            Console.WriteLine(string.Format("SOLUTION: ({0}).", string.Join(", ", solut_best.GetSolut())));
            Console.WriteLine("COST: " + solut_best.GetCost());
        }
    }
}
