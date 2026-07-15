using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace MLP {
    class tSolution {
        private List<int> s;
        private double [] seq;
        private double cost;
        private int d;

        public tSolution(int dimen, double c) {
            seq = new double [(dimen+1)*(dimen+1)*3];
            d = dimen;
            cost = c;
        }

        public void StoreSolut(List<int> sl) {s = sl;}
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double GetSeq(int i, int j, int k) {
            return seq[i*(d+1)*3 + j*3 + k];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void SetSeq(int i, int j, int k, double v) { seq[i*(d+1)*3 + j*3 + k] = v; }

        public double GetCost() {return cost;}
        public void SetCost(double v) {cost = v;}

        public List<int> GetSolut() {return s;}
        public List<int> GetSolutCpy() {return new List<int>(s);}

        public int GetPos(int i) {
            return s[i];
        }

    }

}
