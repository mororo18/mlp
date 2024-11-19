using System.Collections.Generic;
namespace MLP {
    class tSolution {
        private List<int> s;
        private double [][][] seq;
        private double cost;


        public tSolution(int dimen, double c) {
            seq = new double [dimen+1][][];
            for (int i = 0; i < dimen+1; i++) {
                seq[i] = new double [dimen+1][];
                for (int j = 0; j < dimen+1; j++) {
                    seq[i][j] = new double [3];
                }
            }

            cost = c;
        }

        public void StoreSolut(List<int> sl) {s = sl;}
        public double GetSeq(int i, int j, int k) {return seq[i][j][k];}
        public void SetSeq(int i, int j, int k, double v) {seq[i][j][k] = v;}


        public double GetCost() {return cost;}
        public void SetCost(double v) {cost = v;}

        public List<int> GetSolut() {return s;}
        public List<int> GetSolutCpy() {return new List<int>(s);}

    }

}
