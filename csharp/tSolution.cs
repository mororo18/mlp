using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace MLP {
    class tSolution {
        private List<int> s;
        private double [] seq;
        private double cost;
        private int d;

        // `seq` e pinado uma unica vez na construcao (GCHandle, nao um
        // `fixed` por chamada) e o ponteiro fica valido pelo tempo de vida
        // do objeto. solve() cria só 3 instâncias de tSolution, uma vez,
        // reaproveitadas pro programa inteiro — pin-por-chamada em
        // GetSeq/SetSeq foi medido como regressão (1.83x-1.96x mais lento,
        // ver csharp/IMPLEMENTATION.md); pin único no escopo do objeto foi
        // medido como ganho num microbenchmark que aproxima esse cenario.
        // GCHandle nao e liberado explicitamente: soh 3 objetos pinados
        // pelo processo inteiro, que termina logo apos solve().
        private GCHandle seqHandle;
        private unsafe double* seqPtr;

        public unsafe tSolution(int dimen, double c) {
            seq = new double [(dimen+1)*(dimen+1)*3];
            seqHandle = GCHandle.Alloc(seq, GCHandleType.Pinned);
            seqPtr = (double*)seqHandle.AddrOfPinnedObject();
            d = dimen;
            cost = c;
        }

        public void StoreSolut(List<int> sl) {s = sl;}
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe double GetSeq(int i, int j, int k) {
            return seqPtr[i*(d+1)*3 + j*3 + k];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public unsafe void SetSeq(int i, int j, int k, double v) { seqPtr[i*(d+1)*3 + j*3 + k] = v; }

        public double GetCost() {return cost;}
        public void SetCost(double v) {cost = v;}

        public List<int> GetSolut() {return s;}
        public List<int> GetSolutCpy() {return new List<int>(s);}

        public int GetPos(int i) {
            return s[i];
        }

    }

}
