

namespace MLP {
    class tInfo {
        private static double [][] c;
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

        public tInfo(int dimen, double [][] cost, int [] rnd_arr) {
            rnd = rnd_arr;
            rnd_index = 0;
            c = cost;
            dimension = dimen;
        }

        public int GetDimen() {return dimension;}
        public double GetCost(int i, int j) {return c[i][j];}
        public int GetRndCrnt() {return rnd[rnd_index++];}
        
    }

    
}
