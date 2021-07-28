using System;
using System.Collections.Generic;

namespace MLP {
    class GILS_RVND {
        Random rand;
        private double [,] c;
        private double [,,] subseq;

        private int dimension;

        private const int W = 0;
        private const int C = 1;
        private const int T = 2;

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
                cList.Sort((a, b)=> (c[r, a] > c[r, b] ? 1 : -1));
                int range = (int)(((double)cList.Count) * alpha) +1;
                int cN = cList[rand.Next(range)];
                s.Add(cN);
                cList.Remove(cN);
                r = cN;
            }
            s.Add(0);

            return s;
        }

        public void solve(){
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
        }
    }
}
