using System;
using System.Collections.Generic;

namespace MLP {
    class GILS_RVND {
        Random rand;
        private const double EPSILON = 1e-15;
        private double [,] c;
        private double [,,] subseq;

        private int dimension;

        private const int W = 0;
        private const int C = 1;
        private const int T = 2;

        private const int SWAP          = 0;
        private const int REINSERTION   = 1;
        private const int OR_OPT2       = 2;
        private const int OR_OPT3       = 3;
        private const int TWO_OPT       = 4;

        private bool improv_flag;

        private  int                    Iils;
        private const int               Imax = 10;
        private double []  R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
        private const int               R_size = 26;

        public GILS_RVND(){
            Data data = new Data();
            data.loadData();

            dimension = data.getDimension();
            c = new double [dimension, dimension];

            for(int i = 0; i < dimension; i++){
                for(int j = i; j < dimension; j++){
                    c[i, j] = data.getDistance(i, j);
                    c[j, i] = data.getDistance(i, j);
                }
            }

            Iils = (dimension < 100 ? dimension : 100);

            subseq = new double [dimension+1, dimension+1, 3];
            rand = new Random();

        }

        private void subseq_load(List<int> s, double [,,] seq){
            for(int i = 0; i < dimension+1; i++){
                int k = 1 - i - (i != 0 ? 0 : 1);

                seq[i, i, T] = 0.0;
                seq[i, i, C] = 0.0;
                seq[i, i, W] = (i != 0 ? 1.0 : 0.0);

                for(int j = i+1; j < dimension+1; j++){
                    int j_prev = j-1;

                    seq[i, j, T] = c[s[j_prev], s[j]] + seq[i, j_prev, T];
                    seq[i, j, C] = seq[i, j, T] + seq[i, j_prev, C];
                    seq[i, j, W] = j + k;

                    seq[j, i, T] = seq[i, j, T];
                    seq[j, i, C] = seq[i, j, C];
                    seq[j, i, W] = seq[i, j, W];

                    Console.Write(seq[i, j, C] + " ");
                }
                Console.WriteLine();
            }
        }

        private List<int> construction(double alpha){
            var s = new List<int> {0};

            var cList = new List<int>();
            for(int i = 1; i < dimension; i++){
                cList.Add(i);
            }

            int r = 0;
            while(cList.Count > 0){
                cList.Sort((i, j) => (c[r, i] > c[r, j] ? 1 : -1));
                int range = (int)(((double)cList.Count) * alpha) +1;
                int cN = cList[rand.Next(range)];
                s.Add(cN);
                cList.Remove(cN);
                r = cN;
            }
            s.Add(0);

            return s;
        }

        private void swap(List<int> s, int i, int j){
        }

        private void reverse(List<int> s, int i, int j){
        }

        private void reinsert(List<int> s, int i, int j, int pos){
        }

        private void search_swap(List<int> s, double [,,] seq){
        }

        private void search_two_opt(List<int> s, double [,,] seq){
        }

        private void search_reinsertion(List<int> s, double [,,] seq, int opt){
        }

        private void RVND(List<int> s, double [,,] subseq){
            List<int> neighbd_list = new List<int> {SWAP, REINSERTION, OR_OPT2, OR_OPT3, TWO_OPT};

            while(neighbd_list.Count != 0){
                int i_rand = rand.Next(neighbd_list.Count);
                int neighbd = neighbd_list[i_rand];

                improv_flag = false;

                switch(neighbd){
                    case REINSERTION:
                        search_reinsertion(s, subseq, REINSERTION);
                        break;
                    case OR_OPT2:
                        search_reinsertion(s, subseq, OR_OPT2);
                        break;
                    case OR_OPT3:
                        search_reinsertion(s, subseq, OR_OPT3);
                        break;
                    case SWAP:
                        search_swap(s, subseq);
                        break;
                    case TWO_OPT:
                        search_two_opt(s, subseq);
                        break;
                }

                if(improv_flag){
                    neighbd_list.Clear();
                    neighbd_list = new List<int> {SWAP, REINSERTION, OR_OPT2, OR_OPT3, TWO_OPT};
                }else
                    neighbd_list.Remove(i_rand);
            }
        }

        private List<int> perturb(List<int> sl){
            var s = new List<int>(sl);
            return s;
        }


        public void solve(){
            /*
            Console.WriteLine(T);
            var s = new List<int>();
            for(int i = 0; i < dimension; i++){
                s.Add(i);
            }
            s.Add(0);
            s.ForEach(Console.WriteLine);
            subseq_load(s, subseq);
            var opa = construction(0.05);
            Console.WriteLine(string.Format("Inicial: ({0}).", string.Join(", ", opa)));
            */

            double cost_best = Double.MaxValue;
            var s_best = new List<int>();

            for(int i = 0; i < Imax; i++){
                double alpha = R[rand.Next(R_size)];
                var s = construction(alpha);
                var sl = new List<int>(s);

                subseq_load(s, subseq);

                double rvnd_cost_best = subseq[0, dimension, C] - EPSILON;
                double rvnd_cost_crnt;

                Environment.Exit(0);
                int iterILS = 0;
                while(iterILS < Iils){
                    RVND(s, subseq);
                    rvnd_cost_crnt = subseq[0, dimension, C] - EPSILON;
                    if(rvnd_cost_crnt < rvnd_cost_best){
                        rvnd_cost_best = rvnd_cost_crnt;
                        sl.Clear();
                        sl = new List<int>(s);
                        iterILS = 0;
                    }

                    s = perturb(sl);
                    subseq_load(s, subseq);
                    iterILS++;
                }
                subseq_load(sl, subseq);
                double sl_cost = subseq[0, dimension, C] - EPSILON;

                if(sl_cost < cost_best){
                    cost_best = sl_cost;
                    s_best = new List<int>(sl);
                }
            }
        }
    }
}
