
using System.Runtime.CompilerServices;

namespace MLP {
    class tInfo {
        private static double [] c;
        //private static double [][] c;
        private int [] rnd;
        private int rnd_index;

        private int dimension;

        public const int T = 0;
        public const int C = 1;
        public const int W = 2;

        public const int SWAP          = 0;
        public const int REINSERTION   = 1;
        public const int OR_OPT2       = 2;
        public const int OR_OPT3       = 3;
        public const int TWO_OPT       = 4;

        public const double EPSILON = 1e-15;

        private int to_1D(int i, int j) {
            return dimension * i + j;
        }

        public tInfo(int dimen, double [][] cost, int [] rnd_arr) {
            rnd = rnd_arr;
            rnd_index = 0;
            dimension = dimen;
            //c = cost;
            c = new double [dimen*dimen];
            for (int i = 0; i < dimen; i++) {
                for (int j = i; j < dimen; j++) {
                    c[to_1D(i,j)]  = cost[i][j];
                    c[to_1D(j,i)]  = cost[i][j];
                }
                
            }
        }

        public int GetDimen() {return dimension;}
        public double GetCost(int i, int j) {

            //var f = ref c[0];
            double ret = Unsafe.Add(ref c[0], to_1D(i, j));
            return ret;
        }
        //public double GetCost(int i, int j) {return c[i][j];}
        public int GetRndCrnt() {return rnd[rnd_index++];}
        
    }

    
}
