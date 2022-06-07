using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;

namespace MLP {
    class GILS_RVND {
        Random rand;
        private const double EPSILON = 1e-15;
        private static double [][] c;
        //private double [,] c;
        private static double [][][] subseq;
        //private double [,,] subseq;
        private int [] rnd;
        private int rnd_index;

        private int dimension;

        private const int T = 0;
        private const int C = 1;
        private const int W = 2;

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

        // exec time vars
        private long t_construction = 0;
        private long t_swap = 0;
        private long t_two_opt = 0;
        private long t_or_opt2 = 0;
        private long t_or_opt3 = 0;
        private long t_reinsertion = 0;
        private long t_subseq = 0;

        public GILS_RVND(){
            Data data = new Data();
            data.loadData();

            dimension = data.getDimension();
            c = new double [dimension][];
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

            //subseq = new double [dimension+1, dimension+1, 3];
            subseq = new double [dimension+1][][];
            for(int i = 0; i < dimension+1; i++){
                subseq[i] = new double [dimension+1][];

                for(int j = 0; j < dimension+1; j++){
                    subseq[i][j] = new double [3];
                }
            }

            rnd = data.GetRnd();
            rnd_index = 0;

            rand = new Random();

        }

        private void subseq_load(List<int> s, double [][][] seq){
            //t_subseq -= Stopwatch.GetTimestamp();
            for(int i = 0; i < dimension+1; i++){
                int k = 1 - i - (i != 0 ? 0 : 1);

                seq[i][ i][ T] = 0.0;
                seq[i][ i][ C] = 0.0;
                seq[i][ i][ W] = (i != 0 ? 1.0 : 0.0);

                for(int j = i+1; j < dimension+1; j++){
                    int j_prev = j-1;

                    seq[i][ j][ T] = c[s[j_prev]][ s[j]] + seq[i][ j_prev][ T];
                    seq[i][ j][ C] = seq[i][ j][ T] + seq[i][ j_prev][ C];
                    seq[i][ j][ W] = j + k;

                    // a bit faster without it
                    //seq[j, i, T] = seq[i, j, T];
                    //seq[j, i, C] = seq[i, j, C];
                    //seq[j, i, W] = seq[i, j, W];

                    //Console.Write(seq[i, j, C] + " ");
                }
                //Console.WriteLine();
            }
            //t_subseq += Stopwatch.GetTimestamp();
        }

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

        private List<int> construction(double alpha){
            var s = new List<int> {0};

            var cList = new List<int>();
            for(int i = 1; i < dimension; i++){
                cList.Add(i);
            }

            int r = 0;
            while(cList.Count > 0){
                // bug ae (geralmente apos muitas comparacoes)
                //cList.Sort((int i, int j) => (c[r, i]).Equals(c[r, j]));
                cList = cList.OrderBy(i => c[r][ i]).ToList();
                //sort(cList, r);
                //cList.Sort((int i, int j) => (c[r, i] > c[r, j] ? 1 : -1));
                int range = (int)(((double)cList.Count) * alpha) +1;

                int r_value = rand.Next(range);
                r_value = rnd[rnd_index++];

                int cN = cList[r_value];
                s.Add(cN);
                cList.Remove(cN);
                r = cN;
            }
            s.Add(0);

            return s;
        }

        private void swap(List<int> s, int i, int j){
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
            if(i < pos){
                s.InsertRange(pos, s.GetRange(i, sz));
                s.RemoveRange(i, sz);
            }else{
                s.InsertRange(pos, s.GetRange(i, sz));
                s.RemoveRange(i+sz, sz);
            }
            //Console.WriteLine(string.Format("Reinsert: ({0}).", string.Join(", ", s)));
        }

        private void search_swap(List<int> s, double [][][] seq){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1, cost_concat_2,
                   cost_concat_3, cost_concat_4;
            int I = -1;
            int J = -1;

            for(int i = 1; i < dimension-1; i++){
                int i_prev = i - 1;
                int i_next = i + 1;

                cost_concat_1 = seq[0][ i_prev][ T] + c[s[i_prev]][ s[i_next]];
                cost_concat_2 = cost_concat_1 + seq[i][ i_next][ T] + c[s[i]][ s[i_next+1]];

                cost_new = seq[0][ i_prev][ C]
                        + seq[i][ i_next][ W]             * (cost_concat_1) + c[s[i_next]][ s[i]]
                        + seq[i_next+1][ dimension][ W]   * (cost_concat_2) + seq[i_next+1][ dimension][ C];

                if(cost_new < cost_best){
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = i_next;
                }

                for(int j = i_next+1; j < dimension; j++){
                    int j_next = j+1;
                    int j_prev = j-1;

                    cost_concat_1 = seq[0][ i_prev][ T] + c[s[i_prev]][ s[j]];
                    cost_concat_2 = cost_concat_1 + c[s[j]][ s[i_next]];
                    cost_concat_3 = cost_concat_2 + seq[i_next][ j_prev][ T] + c[s[j_prev]][ s[i]];
                    cost_concat_4 = cost_concat_3 + c[s[i]][ s[j_next]];

                    cost_new = seq[0][ i_prev][ C]
                            + cost_concat_1
                            + seq[i_next][ j_prev][ W] * cost_concat_2 + seq[i_next][ j_prev][ C]
                            + cost_concat_3
                            + seq[j_next][ dimension][ W] * cost_concat_4 + seq[j_next][ dimension][ C];
                    
                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if(cost_best < seq[0][ dimension][ C] - EPSILON){
                swap(s, I, J);
                //Console.WriteLine("Swap");
                //Console.WriteLine(cost_best);
                subseq_load(s, seq);
                //Console.WriteLine(seq[0, dimension, C]);
                improv_flag = true;
            }
        }

        private void search_two_opt(List<int> s, double [][][] seq){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1,
                   cost_concat_2;
            int I = -1;
            int J = -1;

            for(int i = 1; i < dimension-1; i++){
                int i_prev = i-1;

                double rev_seq_cost = seq[i][ i+1][ T];

                for(int j = i+2; j < dimension; j++){
                    int j_next = j+1;
                    //int j_prev = j-1;

                    rev_seq_cost += c[s[j-1]][ s[j]] * (seq[i][ j][ W]-1.0);

                    cost_concat_1 = seq[0][ i_prev][ T] + c[s[j]][ s[i_prev]];
                    cost_concat_2 = cost_concat_1 + seq[i][ j][ T] + c[s[j_next]][ s[i]];

                    cost_new = seq[0][ i_prev][ C]
                            + seq[i][ j][ W]              * cost_concat_1 + rev_seq_cost
                            + seq[j_next][ dimension][ W] * cost_concat_2 + seq[j_next][ dimension][ C];

                    /*
                    Console.WriteLine("REV_cost " + rev_seq_cost);
                    Console.WriteLine("concat 1 " + cost_concat_1);
                    Console.WriteLine("concat 2 "  + cost_concat_2);
                    Console.WriteLine(cost_new);
                    */
                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i;
                        J = j;
                    }
                    //Environment.Exit(0);
                }
            }

            if(cost_best < seq[0][ dimension][ C] - EPSILON){
                reverse(s, I, J);
                //Console.WriteLine("Reverse " + I + " " + J);
                //Console.WriteLine(cost_best);
                subseq_load(s, seq);
                //Console.WriteLine(seq[0, dimension, C]);
                improv_flag = true;
            }
            //Environment.Exit(0);
        }

        private void search_reinsertion(List<int> s, double [][][] seq, int opt){
            double cost_best = Double.MaxValue;
            double cost_new;
            double cost_concat_1, cost_concat_2,
                   cost_concat_3;
            int I = -1;
            int J = -1;
            int POS = -1;

            for(int i = 1; i < dimension - opt + 1; i++){
                int j = opt + i - 1;
                int i_prev = i - 1;
                int j_next = j + 1;

                for(int k = 0; k < i_prev; k++){
                    int k_next = k+1;

                    cost_concat_1 = seq[0][ k][ T] + c[s[k]][ s[i]];
                    cost_concat_2 = cost_concat_1 + seq[i][ j][ T] + c[s[j]][ s[k_next]];
                    cost_concat_3 = cost_concat_2 + seq[k_next][ i_prev][ T] + c[s[i_prev]][ s[j_next]];

                    cost_new = seq[0][ k][ C]
                            + seq[i][ j][ W]              * cost_concat_1 + seq[i][ j][ C]
                            + seq[k_next][ i_prev][ W]    * cost_concat_2 + seq[k_next][ i_prev][ C]
                            + seq[j_next][ dimension][ W] * cost_concat_3 + seq[j_next][ dimension][ C];

                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }

                for(int k = i+opt; k < dimension ; k++){
                    int k_next = k+1;

                    cost_concat_1 = seq[0][ i_prev][ T] + c[s[i_prev]][ s[j_next]];;
                    cost_concat_2 = cost_concat_1 + seq[j_next][ k][ T] + c[s[k]][ s[i]];
                    cost_concat_3 = cost_concat_2 + seq[i][ j][ T] + c[s[j]][ s[k_next]];

                    cost_new = seq[0][ i_prev][ C]
                            + seq[j_next][ k][ W]         * cost_concat_1 + seq[j_next][ k][ C]
                            + seq[i][ j][ W]              * cost_concat_2 + seq[i][ j][ C]
                            + seq[k_next][ dimension][ W] * cost_concat_3 + seq[k_next][ dimension][ C];

                    if(cost_new < cost_best){
                        cost_best = cost_new - EPSILON;
                        I = i; 
                        J = j;
                        POS = k;
                    }
                }
            }

            if(cost_best < seq[0][ dimension][ C] - EPSILON){
                //Console.WriteLine("Reinsertion");
                //Console.WriteLine(cost_best);
                reinsert(s, I, J, POS+1);
                subseq_load(s, seq);
                //Console.WriteLine(seq[0, dimension, C]);
                improv_flag = true;
            }

        }

        private void RVND(List<int> s, double [][][] subseq){
            //List<int> neighbd_list = new List<int> {TWO_OPT};
            List<int> neighbd_list = new List<int> {SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3};
            var t = new List<int>();
            for(int i = 0; i < dimension; i++)
                t.Add(i);
            t.Add(0);

            //subseq_load(t, subseq);

            while(neighbd_list.Count != 0){
                int i_rand = rand.Next(neighbd_list.Count);
                i_rand = rnd[rnd_index++];

                int neighbd = neighbd_list[i_rand];

                improv_flag = false;
                //Console.WriteLine(string.Format("T: ({0}).", string.Join(", ", neighbd_list)));

                switch(neighbd){
                    case REINSERTION:
                        //t_reinsertion -= Stopwatch.GetTimestamp();
                        search_reinsertion(s, subseq, REINSERTION);
                        //t_reinsertion += Stopwatch.GetTimestamp();
                        break;
                    case OR_OPT2:
                        //t_or_opt2 -= Stopwatch.GetTimestamp();
                        search_reinsertion(s, subseq, OR_OPT2);
                        //t_or_opt2 += Stopwatch.GetTimestamp();
                        break;
                    case OR_OPT3:
                        //t_or_opt3 -= Stopwatch.GetTimestamp();
                        search_reinsertion(s, subseq, OR_OPT3);
                        //t_or_opt3 += Stopwatch.GetTimestamp();
                        break;
                    case SWAP:
                        //t_swap -= Stopwatch.GetTimestamp();
                        search_swap(s, subseq);
                        //t_swap += Stopwatch.GetTimestamp();
                        break;
                    case TWO_OPT:
                        //t_two_opt -= Stopwatch.GetTimestamp();
                        search_two_opt(s, subseq);
                        //t_two_opt += Stopwatch.GetTimestamp();
                        break;
                }

                //Environment.Exit(0);

                if(improv_flag){
                    //neighbd_list.Clear();
                    neighbd_list = new List<int> {SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3};
                    //Console.WriteLine("not Removed " + i_rand);
                }else{
                    //Console.WriteLine("Removed " + i_rand);
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
                int max = sl.Count - 1 - size_max;
                A_start = rand.Next(max) + 1;
                A_end = A_start + rand.Next(size_max - size_min + 1) + size_min;
                B_start = rand.Next(max) + 1;
                B_end = B_start + rand.Next(size_max - size_min + 1) + size_min;

                A_start = rnd[rnd_index++];
                A_end = A_start + rnd[rnd_index++];
                
                B_start = rnd[rnd_index++];
                B_end = B_start + rnd[rnd_index++];
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
                int index = rand.Next(R_size);
                index = rnd[rnd_index++];
                double alpha = R[index];

                Console.WriteLine("[+] Local Search " + (i+1));
                Console.WriteLine("\t[+] Constructing Inital Solution..");

                //t_construction -= Stopwatch.GetTimestamp();
                var s = construction(alpha);

                /*
                for (int j = 0; j < dimension; j++) {
                    Console.Write(s[j] + " ");
                }
                Console.WriteLine();
                */
                //t_construction += Stopwatch.GetTimestamp();
                var sl = new List<int>(s);

                subseq_load(s, subseq);

                double rvnd_cost_best = subseq[0][ dimension][ C] - EPSILON;
                double rvnd_cost_crnt;

                //Environment.Exit(0);
                Console.WriteLine("\t[+] Looking for the best Neighbor..");
                int iterILS = 0;
                while(iterILS < Iils){
                    RVND(s, subseq);
                    rvnd_cost_crnt = subseq[0][ dimension][ C] - EPSILON;
                    if(rvnd_cost_crnt < rvnd_cost_best){
                        rvnd_cost_best = rvnd_cost_crnt;
                        //sl.Clear();
                        sl = new List<int>(s);
                        iterILS = 0;
                    }

                    s = perturb(sl);
                    subseq_load(s, subseq);
                    iterILS++;
                }
                subseq_load(sl, subseq);
                double sl_cost = subseq[0][ dimension][ C] - EPSILON;

                if(sl_cost < cost_best){
                    cost_best = sl_cost;
                    s_best = new List<int>(sl);
                }
                Console.WriteLine("\tCurrent best solution cost: "+cost_best);
            }
            Console.WriteLine(string.Format("SOLUTION: ({0}).", string.Join(", ", s_best)));
            Console.WriteLine("COST: " + cost_best);
            Console.WriteLine("Construction: " + t_construction/10e6);
            Console.WriteLine("Swap: " + t_swap/10e6);
            Console.WriteLine("2_opt: " + t_two_opt/10e6);
            Console.WriteLine("or_opt2: " + t_or_opt2/10e6);
            Console.WriteLine("or_opt3: " + t_or_opt3/10e6);
            Console.WriteLine("reinseriton: " + t_reinsertion/10e6);
            Console.WriteLine("subseq_load: " + t_subseq/10e6);
        }
    }
}
