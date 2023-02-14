using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace MLP {
    class tSolution {
        private List<int> s;
        private double [] seq;
        //private double [][] seq;
        private double cost;
        private int d;

        /*public tSolution(int dimen, double c) {
            seq = new double [dimen+1][][];
            for (int i = 0; i < dimen+1; i++) {
                seq[i] = new double [dimen+1][];
                for (int j = 0; j < dimen+1; j++) {
                    seq[i][j] = new double [3];
                }
            }

            cost = c;
        }*/

        public tSolution(int dimen, double c) {
            
            /*
            seq = new double [dimen+1][];
            for (int i = 0; i < dimen+1; i++) {
                seq[i] = new double [(dimen+1)*3];
            }
            */

            seq = new double[(dimen + 1)*(dimen+1)*3];
            d = dimen;

            cost = c;
        }

        public void StoreSolut(List<int> sl) {s = sl;}
        public double GetSeq(int i, int j, int k) {
            /*
            ref var ptr = ref seq[0];
            var nptr = Unsafe.Add(ref ptr, i);
            ref var nn = ref nptr;
            var ret  = Unsafe.Add(ref nn, 3*j + k);
            */

            //ref var ptr = ref seq[0];
            double ret = Unsafe.Add(ref seq[0], 3*((d+1)*i + j) + k);
            return ret;

            //return seq[i*3*(d+1)+j*3 + k];
            //return seq[i][j*3 + k];
        }

        public void SetSeq(int i, int j, int k, double v) {Unsafe.Add(ref seq[0], 3*((d+1)*i + j) + k) = v;}
        //public void SetSeq(int i, int j, int k, double v) {seq[i*3*(d+1)+j*3 + k] = v;}

        //public void SetSeq(int i, int j, int k, double v) { seq[i][j * 3 + k] = v; }


        /*        
        public double GetSeq(int i, int j, int k) {return seq[i][j][k];}
        public void SetSeq(int i, int j, int k, double v) {seq[i][j][k] = v;}
        */

        public double GetCost() {return cost;}
        public void SetCost(double v) {cost = v;}

        public List<int> GetSolut() {return s;}
        public List<int> GetSolutCpy() {return new List<int>(s);}

        public int GetPos(int i) {
            return Unsafe.Add(ref s[0], i);
        }

    }

}
