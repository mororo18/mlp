using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;

namespace MLP {
    class GILS_RVND {
        tInfo info;

        private  int                    Iils;
        private const int               Imax = 10;
        private double []  R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
        private const int               R_size = 26;

        public GILS_RVND(){
            Data data = new Data();
            data.loadData();

            Console.WriteLine(string.Format("dimesion: {0}.", data.getDimension()));

            int dimension = data.getDimension();
            double [][] c = new double [dimension][];
            for(int i = 0; i < dimension; i++){
                c[i] = new double [dimension];
            }

            for(int i = 0; i < dimension; i++){
                for(int j = i; j < dimension; j++){
                    c[i][j] = data.getDistance(i, j);
                    c[j][i] = data.getDistance(j, i);
                }
            }

            Iils = (dimension < 100 ? dimension : 100);


            int [] rnd = data.GetRnd();

            info = new tInfo(dimension, c, rnd);
        }

        private void subseq_load(tSolution solut, tInfo info){
            var solutSolut = solut.GetSolut();  
            for(int i = 0; i < info.GetDimen()+1; i++){
                int k = 1 - i - (i != 0 ? 0 : 1);

                solut.SetSeq(i, i, tInfo.T, 0.0);
                solut.SetSeq(i, i, tInfo.C, 0.0);
                solut.SetSeq(i, i, tInfo.W, (i != 0 ? 1.0 : 0.0));

                for(int j = i+1; j < info.GetDimen()+1; j++){
                    int j_prev = j-1;

                    solut.SetSeq(i, j, tInfo.T, info.GetCost(solutSolut[j_prev], solutSolut[j]) + solut.GetSeq(i, j_prev, tInfo.T));
                    solut.SetSeq(i, j, tInfo.C, solut.GetSeq(i, j, tInfo.T) + solut.GetSeq(i, j_prev, tInfo.C));
                    solut.SetSeq(i, j, tInfo.W, j + k);

                }
            }

            solut.SetCost(solut.GetSeq(0, info.GetDimen(), tInfo.C));
        }

        /*
        private void sort(List<int> arr, int r) {
            for (int i = 0; i < arr.Count; i++) {
                for (int j = 0; j < arr.Count - i-1; j++) {
                    //Console.WriteLine(arr[j+1]);
                    if (c[r][arr[j]] > c[r][arr[j+1]]) {
                        int tmp = arr[j];
                        arr[j] = arr[j+1];
                        arr[j+1] = tmp;
                    }
                }
            }
        }
        */

        private void sort(List<int> arr, int r, tInfo info)
        {
            Quicksort(arr, 0, arr.Count - 1, r, info);
        }

        private void Quicksort(List<int> arr, int left, int right, int r, tInfo info)
        {
            if (left < right)
            {
                int pivot = Partition(arr, left, right, r, info);
                Quicksort(arr, left, pivot - 1, r, info);
                Quicksort(arr, pivot + 1, right, r, info);
            }
        }

        private int Partition(List<int> arr, int left, int right, int r, tInfo info)
        {
            int pivotIndex = right;
            double pivotValue = info.GetCost(r, arr[pivotIndex]);
            int i = left - 1;

            for (int j = left; j < right; j++)
            {
                if (info.GetCost(r, arr[j]) < pivotValue)
                {
                    i++;
                    (arr[i], arr[j]) = (arr[j], arr[i]); // Troca os elementos usando tuple swap
                }
            }
            (arr[i + 1], arr[right]) = (arr[right], arr[i + 1]); // Troca o pivô
            return i + 1;
        }

        

        private List<int> construction(double alpha, tInfo info){
            var s = new List<int> {0};

            var cList = new List<int>();
            for(int i = 1; i < info.GetDimen(); i++){
                cList.Add(i);
            }

            int r = 0;
            while(cList.Count > 0){

                sort(cList, r, info);

                int r_value = info.GetRndCrnt();

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

        private bool search_swap(tSolution solut, tInfo info){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1, cost_concat_2,
                   cost_concat_3, cost_concat_4;
            int I = -1;
            int J = -1;

            var solutSolut = solut.GetSolut();
            var dimen = info.GetDimen();

            for(int i = 1; i < dimen-1; i++){
                int i_prev = i - 1;
                int i_next = i + 1;

                cost_concat_1 = solut.GetSeq(0, i_prev, tInfo.T) + info.GetCost(solut.GetPos(i_prev),  solut.GetPos(i_next));
                cost_concat_2 = cost_concat_1 + solut.GetSeq(i, i_next, tInfo.T) + info.GetCost(solutSolut[i], solutSolut[i_next+1]);

                cost_new = solut.GetSeq(0, i_prev, tInfo.C)
                        + solut.GetSeq(i, i_next, tInfo.W)             * (cost_concat_1) + info.GetCost(solutSolut[i_next], solutSolut[i])
                        + solut.GetSeq(i_next+1, dimen, tInfo.W)   * (cost_concat_2) + solut.GetSeq(i_next+1, dimen, tInfo.C);

                if(cost_new < cost_best){
                    cost_best = cost_new - tInfo.EPSILON;
                    I = i;
                    J = i_next;
                }

                for(int j = i_next+1; j < dimen; j++){
                    int j_next = j+1;
                    int j_prev = j-1;

                    cost_concat_1 = solut.GetSeq(0, i_prev, tInfo.T) + info.GetCost(solutSolut[i_prev], solutSolut[j]);
                    cost_concat_2 = cost_concat_1 + info.GetCost(solutSolut[j], solutSolut[i_next]);
                    cost_concat_3 = cost_concat_2 + solut.GetSeq(i_next, j_prev, tInfo.T) + info.GetCost(solutSolut[j_prev], solutSolut[i]);
                    cost_concat_4 = cost_concat_3 + info.GetCost(solutSolut[i], solutSolut[j_next]);

                    cost_new = solut.GetSeq(0, i_prev, tInfo.C)
                            + cost_concat_1
                            + solut.GetSeq(i_next, j_prev, tInfo.W) * cost_concat_2 + solut.GetSeq(i_next, j_prev, tInfo.C)
                            + cost_concat_3
                            + solut.GetSeq(j_next, dimen, tInfo.W) * cost_concat_4 + solut.GetSeq(j_next, dimen, tInfo.C);
                    
                    if(cost_new < cost_best){
                        cost_best = cost_new - tInfo.EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if(cost_best < solut.GetCost() - tInfo.EPSILON){
              //if (I == -1 && J == -1) {
              //    Console.WriteLine(string.Format("Swap: ({0}, {1}, {2}).", cost_best, solut.GetCost(), -cost_best+solut.GetCost()));
              //}
                swap(solutSolut, I, J);
                subseq_load(solut, info);
                return true;
            }

            return false;
        }

        private bool search_two_opt(tSolution solut, tInfo info){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1,
                   cost_concat_2;
            int I = -1;
            int J = -1;
            var solutSolut = solut.GetSolut();

            var dimen = info.GetDimen();

            for (int i = 1; i < dimen-1; i++){
                int i_prev = i-1;

                double rev_seq_cost = solut.GetSeq(i, i+1, tInfo.T);

                for(int j = i+2; j < dimen; j++){
                    int j_next = j+1;
                    //int j_prev = j-1;

                    rev_seq_cost += info.GetCost(solutSolut[j-1], solutSolut[j]) * (solut.GetSeq(i,  j, tInfo.W)-1.0);

                    cost_concat_1 = solut.GetSeq(0, i_prev, tInfo.T) + info.GetCost(solutSolut[j], solutSolut[i_prev]);
                    cost_concat_2 = cost_concat_1 + solut.GetSeq(i, j,  tInfo.T) + info.GetCost(solutSolut[j_next], solutSolut[i]);

                    cost_new = solut.GetSeq(0, i_prev,  tInfo.C)
                            + solut.GetSeq(i, j, tInfo.W)              * cost_concat_1 + rev_seq_cost
                            + solut.GetSeq(j_next, dimen, tInfo.W) * cost_concat_2 + solut.GetSeq(j_next, dimen, tInfo.C);

                    if(cost_new < cost_best){
                        cost_best = cost_new - tInfo.EPSILON;
                        I = i;
                        J = j;
                    }

                    //Environment.Exit(0);
                }
            }

            if(cost_best < solut.GetCost() - tInfo.EPSILON){
                reverse(solutSolut, I, J);
                subseq_load(solut, info);
                return true;
            }

            return false;
        }

        private bool search_reinsertion(tSolution solut, tInfo info, int opt){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1, cost_concat_2,
                   cost_concat_3;
            int I = -1;
            int J = -1;
            int POS = -1;
            var solutSolut = solut.GetSolut();

            var dimen = info.GetDimen();

            for (int i = 1; i < dimen - opt + 1; i++){
                int j = opt + i - 1;
                int i_prev = i - 1;
                int j_next = j + 1;

                for(int k = 0; k < i_prev; k++){
                    int k_next = k+1;

                    cost_concat_1 = solut.GetSeq(0, k, tInfo.T) + info.GetCost(solutSolut[k], solutSolut[i]);
                    cost_concat_2 = cost_concat_1 + solut.GetSeq(i, j, tInfo.T) + info.GetCost(solutSolut[j], solutSolut[k_next]);
                    cost_concat_3 = cost_concat_2 + solut.GetSeq(k_next, i_prev, tInfo.T) + info.GetCost(solutSolut[i_prev], solutSolut[j_next]);

                    cost_new = solut.GetSeq(0, k, tInfo.C)
                            + solut.GetSeq(i, j, tInfo.W)              * cost_concat_1 + solut.GetSeq(i, j, tInfo.C)
                            + solut.GetSeq(k_next, i_prev, tInfo.W)    * cost_concat_2 + solut.GetSeq(k_next, i_prev, tInfo.C)
                            + solut.GetSeq(j_next, dimen, tInfo.W) * cost_concat_3 + solut.GetSeq(j_next, dimen, tInfo.C);

                    if(cost_new < cost_best){
                        cost_best = cost_new - tInfo.EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }

                for(int k = i+opt; k < dimen; k++){
                    int k_next = k+1;

                    cost_concat_1 = solut.GetSeq(0, i_prev, tInfo.T) + info.GetCost(solutSolut[i_prev], solutSolut[j_next]);
                    cost_concat_2 = cost_concat_1 + solut.GetSeq(j_next, k, tInfo.T) + info.GetCost(solutSolut[k], solutSolut[i]);
                    cost_concat_3 = cost_concat_2 + solut.GetSeq(i, j, tInfo.T) + info.GetCost(solutSolut[j], solutSolut[k_next]);

                    cost_new = solut.GetSeq(0, i_prev, tInfo.C)
                            + solut.GetSeq(j_next, k, tInfo.W)         * cost_concat_1 + solut.GetSeq(j_next, k, tInfo.C)
                            + solut.GetSeq(i, j, tInfo.W)              * cost_concat_2 + solut.GetSeq(i, j, tInfo.C)
                            + solut.GetSeq(k_next, dimen, tInfo.W) * cost_concat_3 + solut.GetSeq(k_next, dimen, tInfo.C);

                    if(cost_new < cost_best){
                        cost_best = cost_new - tInfo.EPSILON;
                        I = i; 
                        J = j;
                        POS = k;
                    }
                }
            }

            if(cost_best < solut.GetCost() - tInfo.EPSILON){
                //Console.WriteLine("Reinsertion");
                //Console.WriteLine(cost_best);
                reinsert(solutSolut, I, J, POS+1);
                subseq_load(solut, info);
                //Console.WriteLine(seq[0, dimension, C]);
                return true;
            }

            return false;

        }

        private void RVND(tSolution solut, tInfo info){
            List<int> neighbd_list = new List<int> {tInfo.SWAP, tInfo.TWO_OPT, tInfo.REINSERTION, tInfo.OR_OPT2, tInfo.OR_OPT3};

            /*
            var t = new List<int>();
            for(int i = 0; i < dimension; i++)
                t.Add(i);
            t.Add(0);
            */

            while(neighbd_list.Count != 0){
                int i_rand = info.GetRndCrnt();

                int neighbd = neighbd_list[i_rand];

                bool improve = false;

                switch(neighbd){
                    case tInfo.REINSERTION:
                        improve = search_reinsertion(solut, info, tInfo.REINSERTION);
                        break;
                    case tInfo.OR_OPT2:
                        improve = search_reinsertion(solut, info, tInfo.OR_OPT2);
                        break;
                    case tInfo.OR_OPT3:
                        improve = search_reinsertion(solut, info, tInfo.OR_OPT3);
                        break;
                    case tInfo.SWAP:
                        improve = search_swap(solut, info);
                        break;
                    case tInfo.TWO_OPT:
                        improve = search_two_opt(solut, info);
                        break;
                }

                if(improve){
                    neighbd_list = new List<int> {tInfo.SWAP, tInfo.TWO_OPT, tInfo.REINSERTION, tInfo.OR_OPT2, tInfo.OR_OPT3};
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

                A_start = info.GetRndCrnt();
                A_end = A_start + info.GetRndCrnt();
                
                B_start = info.GetRndCrnt();
                B_end = B_start + info.GetRndCrnt();
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
            var solut_best = new tSolution(info.GetDimen(), Double.MaxValue);
            var solut_crnt = new tSolution(info.GetDimen(), 0.0);
            var solut_partial = new tSolution(info.GetDimen(), 0.0);

            for(int i = 0; i < Imax; i++){
                int index = info.GetRndCrnt();
                double alpha = R[index];

                Console.WriteLine("[+] Local Search " + (i+1));
                Console.WriteLine("\t[+] Constructing Inital Solution..");

                solut_crnt.StoreSolut(construction(alpha, info));

                subseq_load(solut_crnt, info);

                solut_partial.StoreSolut(solut_crnt.GetSolutCpy());
                solut_partial.SetCost(solut_crnt.GetCost());

                Console.WriteLine("Construction: " + solut_crnt.GetCost());

                Console.WriteLine("\t[+] Looking for the best Neighbor..");
                int iterILS = 0;
                while(iterILS < Iils){
                    RVND(solut_crnt, info);
                    if(solut_crnt.GetCost() < solut_partial.GetCost()){
                        solut_partial.StoreSolut(solut_crnt.GetSolutCpy());
                        solut_partial.SetCost(solut_crnt.GetCost());

                        iterILS = 0;
                    }

                    solut_crnt.StoreSolut(perturb(solut_partial.GetSolut()));
                    subseq_load(solut_crnt, info);
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
