import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Random;

public class Main {

    static final int TABLE_SZ = 26;

    static final double DBL_EPSILON = 2.220446049250313e-16;
    static final double DBL_MAX = Double.MAX_VALUE;

    static class tData {
        double[][] cost;
        int dimen;
        int[] rnd;
        int rnd_index;
    }

    // Neighborhood
    static final int REINSERTION = 1;
    static final int OR_OPT_2 = 2;
    static final int OR_OPT_3 = 3;
    static final int SWAP = 4;
    static final int TWO_OPT = 5;

    // the tSeqInfo (T, C, W) matrix is kept flattened, one array per field
    static class tSolution {
        int[] s;
        double[] seqT;
        double[] seqC;
        double[] seqW;
        double cost;
    }

    static int to_1D(int x, int y, int size) {
        int a = x * (size + 1) + y;

        return a;
    }

    static tSolution Solution_init(tData data) {
        tSolution solut = new tSolution();
        solut.s = new int[data.dimen + 1];

        int sz = (data.dimen + 1) * (data.dimen + 1);
        solut.seqT = new double[sz];
        solut.seqC = new double[sz];
        solut.seqW = new double[sz];

        solut.cost = DBL_MAX;

        return solut;
    }

    static void Solution_cpy(tSolution src, tSolution tgt, tData data) {

        tgt.s = src.s.clone();
        tgt.cost = src.cost;

    }

    static final double[] R_TABLE = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                                     0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};

    static double R_table(int i) {
        return R_TABLE[i];
    }

    static void print_s(int[] s) {

        for (int i = 0; i < s.length; i++)
            System.out.print((s[i] + 1) + " ");
        System.out.println();
    }

    static int partition(int[] arr, int left, int right, tData data, int r) {
        int pivot = arr[right];
        int i = left - 1;
        for (int j = left; j < right; j++) {
            if (data.cost[r][arr[j]] < data.cost[r][pivot]) {
                i++;
                int tmp = arr[i];
                arr[i] = arr[j];
                arr[j] = tmp;
            }
        }
        int tmp = arr[i + 1];
        arr[i + 1] = arr[right];
        arr[right] = tmp;
        return i + 1;
    }

    static void quicksort(int[] arr, int left, int right, tData data, int r) {
        if (left < right) {
            int pivot = partition(arr, left, right, data, r);
            quicksort(arr, left, pivot - 1, data, r);
            quicksort(arr, pivot + 1, right, data, r);
        }
    }

    static void sort(int[] arr, int size, int r, tData data) {
        quicksort(arr, 0, size - 1, data, r);
    }

    static int[] construct(final double alpha, tData data) {

        int[] s = new int[data.dimen + 1];
        int s_size = 0;
        s[s_size++] = 0;

        int cL_size = data.dimen - 1;
        int[] cL = new int[cL_size];
        for (int i = 1; i < data.dimen; ++i) {
            cL[i - 1] = i;
        }

        int r = 0;
        while (cL_size > 0) {
            sort(cL, cL_size, r, data);

            int index = data.rnd[data.rnd_index++];

            int c = cL[index];
            s[s_size++] = c;
            r = c;
            System.arraycopy(cL, index + 1, cL, index, cL_size - index - 1);
            cL_size--;
        }

        s[s_size++] = 0;

        return s;
    }

    static void swap(int[] vec, int i, int j) {
        int tmp = vec[i];
        vec[i] = vec[j];
        vec[j] = tmp;
    }

    static void reverse(int[] vec, int i, int j) {
        while (i < j) {
            int tmp = vec[i];
            vec[i] = vec[j];
            vec[j] = tmp;
            i++;
            j--;
        }
    }

    static void reinsert(int[] vec, int i, int j, int pos) {
        int len = j - i + 1;
        int[] seq = new int[len];
        System.arraycopy(vec, i, seq, 0, len);
        if (pos < i) {
            System.arraycopy(vec, pos, vec, pos + len, i - pos);
            System.arraycopy(seq, 0, vec, pos, len);
        } else {
            System.arraycopy(vec, j + 1, vec, i, pos - j - 1);
            System.arraycopy(seq, 0, vec, pos - len, len);
        }

    }

    static void update_subseq_info_matrix(tSolution solut, tData data) {
        update_subseq_info_matrix(solut, data, 0);
    }

    static void update_subseq_info_matrix(tSolution solut, tData data, int index) {
        int i, j, j_prev, k;
        int from = index;
        boolean t;
        for (i = 0; i < data.dimen + 1; i++) {
            k = 1 - i - (i == 0 ? 1 : 0);
            t = i == from;

            solut.seqT[to_1D(i, i, data.dimen)] = 0.0;
            solut.seqC[to_1D(i, i, data.dimen)] = 0.0;
            solut.seqW[to_1D(i, i, data.dimen)] = (i == 0 ? 0.0 : 1.0);

            for (j = i + 1; j < data.dimen + 1; j++) {
                j_prev = j - 1;

                solut.seqT[to_1D(i, j, data.dimen)] = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seqT[to_1D(i, j_prev, data.dimen)];
                solut.seqC[to_1D(i, j, data.dimen)] = solut.seqT[to_1D(i, j, data.dimen)] + solut.seqC[to_1D(i, j_prev, data.dimen)];
                solut.seqW[to_1D(i, j, data.dimen)] = j + k;
            }
            from += t ? 1 : 0;
        }

        solut.cost = solut.seqC[to_1D(0, data.dimen, data.dimen)];
    }

    static boolean search_swap(tSolution solut, tData data) {
        double cost_new,
            cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
        double cost_best = DBL_MAX;
        int i, j, j_prev, j_next, i_prev, i_next;
        int I = 0;
        int J = 0;

        for (i = 1; i < data.dimen - 1; ++i) {
            i_prev = i - 1;
            i_next = i + 1;

            //consecutive nodes
            cost_concat_1 =                 solut.seqT[to_1D(0, i_prev, data.dimen)] + data.cost[solut.s[i_prev]][solut.s[i_next]];
            cost_concat_2 = cost_concat_1 + solut.seqT[to_1D(i, i_next, data.dimen)]  + data.cost[solut.s[i]][solut.s[i_next + 1]];

            cost_new = solut.seqC[to_1D(0, i_prev, data.dimen)]                                                    +           //       1st subseq
            solut.seqW[to_1D(i, i_next, data.dimen)]               * (cost_concat_1) + data.cost[solut.s[i_next]][solut.s[i]]  +           // concat 2nd subseq
            solut.seqW[to_1D(i_next + 1, data.dimen, data.dimen)]   * (cost_concat_2) + solut.seqC[to_1D(i_next + 1, data.dimen, data.dimen)];   // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = i_next;
            }

            for (j = i_next + 1; j < data.dimen; ++j) {
                j_next = j + 1;
                j_prev = j - 1;

                cost_concat_1 =                 solut.seqT[to_1D(0, i_prev, data.dimen)]       + data.cost[solut.s[i_prev]][solut.s[j]];
                cost_concat_2 = cost_concat_1                           + data.cost[solut.s[j]][solut.s[i_next]];
                cost_concat_3 = cost_concat_2 + solut.seqT[to_1D(i_next, j_prev, data.dimen)]  + data.cost[solut.s[j_prev]][solut.s[i]];
                cost_concat_4 = cost_concat_3                           + data.cost[solut.s[i]][solut.s[j_next]];

                cost_new = solut.seqC[to_1D(0, i_prev, data.dimen)]                                                 +      // 1st subseq
                cost_concat_1 +                                                             // concat 2nd subseq (single node)
                solut.seqW[to_1D(i_next, j_prev, data.dimen)]      * cost_concat_2 + solut.seqC[to_1D(i_next, j_prev, data.dimen)] +      // concat 3rd subseq
                cost_concat_3 +                                                             // concat 4th subseq (single node)
                solut.seqW[to_1D(j_next, data.dimen, data.dimen)] * cost_concat_4 + solut.seqC[to_1D(j_next, data.dimen, data.dimen)];   // concat 5th subseq

                if (cost_new < cost_best) {
                    cost_best = cost_new - DBL_EPSILON;
                    I = i;
                    J = j;
                }
            }
        }

        if (cost_best < solut.cost - DBL_EPSILON) {
            swap(solut.s, I, J);
            update_subseq_info_matrix(solut, data, I);
            return true;
        }

        return false;
    }

    static boolean search_two_opt(tSolution solut, tData data) {
        double cost_new,
            cost_concat_1, cost_concat_2;
        double cost_best = DBL_MAX; // cost_l1, cost_l2;
        double rev_seq_cost;
        int i, j, i_prev, j_next;
        int I = 0;
        int J = 0;

        for (i = 1; i < data.dimen - 1; ++i) {
            i_prev = i - 1;

            rev_seq_cost = solut.seqT[to_1D(i, i + 1, data.dimen)];
            for (j = i + 2; j < data.dimen; ++j) {
                j_next = j + 1;

                rev_seq_cost += data.cost[solut.s[j - 1]][solut.s[j]] * (solut.seqW[to_1D(i, j, data.dimen)] - 1.0);

                cost_concat_1 =                 solut.seqT[to_1D(0, i_prev, data.dimen)]   + data.cost[solut.s[j]][solut.s[i_prev]];
                cost_concat_2 = cost_concat_1 + solut.seqT[to_1D(i, j, data.dimen)]        + data.cost[solut.s[j_next]][solut.s[i]];

                cost_new = solut.seqC[to_1D(0, i_prev, data.dimen)]                                                        +   //  1st subseq
                        solut.seqW[to_1D(i, j, data.dimen)]                * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
                        solut.seqW[to_1D(j_next, data.dimen, data.dimen)] * cost_concat_2 + solut.seqC[to_1D(j_next, data.dimen, data.dimen)];      // concat 3rd subseq

                if (cost_new < cost_best) {
                    cost_best = cost_new - DBL_EPSILON;
                    I = i;
                    J = j;
                }
            }
        }

        if (cost_best < solut.cost - DBL_EPSILON) {
            reverse(solut.s, I, J);
            update_subseq_info_matrix(solut, data);
            return true;
        }

        return false;
    }

    static boolean search_reinsertion(tSolution solut, tData data, final int opt) {
        double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
        double cost_best = DBL_MAX; //, cost_l1, cost_l2, cost_l3;
        int i, j, k, k_next, i_prev, j_next;
        int I = 0;
        int J = 0;
        int POS = 0;

        for (i = 1, j = opt + i - 1; i < data.dimen - opt + 1; ++i, ++j) {
            j_next = j + 1;
            i_prev = i - 1;

            //k -> edges
            for (k = 0; k < i_prev; ++k) {
                k_next = k + 1;

                cost_concat_1 =                 solut.seqT[to_1D(0, k, data.dimen)]            + data.cost[solut.s[k]][solut.s[i]];
                cost_concat_2 = cost_concat_1 + solut.seqT[to_1D(i, j, data.dimen)]            + data.cost[solut.s[j]][solut.s[k_next]];
                cost_concat_3 = cost_concat_2 + solut.seqT[to_1D(k_next, i_prev, data.dimen)]  + data.cost[solut.s[i_prev]][solut.s[j_next]];

                cost_new = solut.seqC[to_1D(0, k, data.dimen)]                                                                   +   //       1st subseq
                    solut.seqW[to_1D(i, j, data.dimen)]               * cost_concat_1 + solut.seqC[to_1D(i, j, data.dimen)]                  +   //  concat 2nd subseq (reinserted seq)
                    solut.seqW[to_1D(k_next, i_prev, data.dimen)]     * cost_concat_2 + solut.seqC[to_1D(k_next, i_prev, data.dimen)]        +   //  concat 3rd subseq
                    solut.seqW[to_1D(j_next, data.dimen, data.dimen)] * cost_concat_3 + solut.seqC[to_1D(j_next, data.dimen, data.dimen)];       // concat 4th subseq

                if (cost_new < cost_best) {
                    cost_best = cost_new - DBL_EPSILON;
                    I = i;
                    J = j;
                    POS = k;
                }
            }

            for (k = i + opt; k < data.dimen; ++k) {
                k_next = k + 1;

                cost_concat_1 =                 solut.seqT[to_1D(0, i_prev, data.dimen)]  + data.cost[solut.s[i_prev]][solut.s[j_next]];
                cost_concat_2 = cost_concat_1 + solut.seqT[to_1D(j_next, k, data.dimen)]  + data.cost[solut.s[k]][solut.s[i]];
                cost_concat_3 = cost_concat_2 + solut.seqT[to_1D(i, j, data.dimen)]       + data.cost[solut.s[j]][solut.s[k_next]];

                cost_new = solut.seqC[to_1D(0, i_prev, data.dimen)]                                                                  +   //       1st subseq
                    solut.seqW[to_1D(j_next, k, data.dimen)]          * cost_concat_1 + solut.seqC[to_1D(j_next, k, data.dimen)]             +   // concat 2nd subseq
                    solut.seqW[to_1D(i, j, data.dimen)]               * cost_concat_2 + solut.seqC[to_1D(i, j, data.dimen)]                  +   // concat 3rd subseq (reinserted seq)
                    solut.seqW[to_1D(k_next, data.dimen, data.dimen)] * cost_concat_3 + solut.seqC[to_1D(k_next, data.dimen, data.dimen)];       // concat 4th subseq


                if (cost_new < cost_best) {
                    cost_best = cost_new - DBL_EPSILON;
                    I = i;
                    J = j;
                    POS = k;
                }
            }
        }

        if (cost_best < solut.cost - DBL_EPSILON) {
            reinsert(solut.s, I, J, POS + 1);
            update_subseq_info_matrix(solut, data, I < POS + 1 ? I : POS + 1);
            return true;
        }

        return false;
    }

    static void RVND(tSolution solut, tData data) {

        int[] neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
        int neighbd_list_size = 5;
        int index;
        int neighbd;
        boolean improve_flag;

        while (neighbd_list_size > 0) {

            index = data.rnd[data.rnd_index++];
            neighbd = neighbd_list[index];
            improve_flag = false;

            switch (neighbd) {
                case REINSERTION:
                    improve_flag = search_reinsertion(solut, data, REINSERTION);
                    break;
                case OR_OPT_2:
                    improve_flag = search_reinsertion(solut, data, OR_OPT_2);
                    break;
                case OR_OPT_3:
                    improve_flag = search_reinsertion(solut, data, OR_OPT_3);
                    break;
                case SWAP:
                    improve_flag = search_swap(solut, data);
                    break;
                case TWO_OPT:
                    improve_flag = search_two_opt(solut, data);
                    break;
            }
            if (improve_flag) {
                neighbd_list = new int[]{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
                neighbd_list_size = 5;
            } else {
                System.arraycopy(neighbd_list, index + 1, neighbd_list, index, neighbd_list_size - index - 1);
                neighbd_list_size--;
            }

        }

    }

    static int[] perturb(tSolution solut, tData data) {
        int[] s = solut.s.clone();
        int A_start = 1;
        int A_end = 1;
        int B_start = 1;
        int B_end = 1;

        int size_max = (int) Math.floor((data.dimen + 1) / 10.0);
        size_max = size_max >= 2 ? size_max : 2;
        int size_min = 2;
        while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

            A_start = data.rnd[data.rnd_index++];
            A_end = A_start + data.rnd[data.rnd_index++];

            B_start = data.rnd[data.rnd_index++];
            B_end = B_start + data.rnd[data.rnd_index++];
        }

        if (A_start < B_start) {
            reinsert(s, B_start, B_end - 1, A_end);
            reinsert(s, A_start, A_end - 1, B_end);
        } else {
            reinsert(s, A_start, A_end - 1, B_end);
            reinsert(s, B_start, B_end - 1, A_end);
        }

        return s;
    }


    static void GILS_RVND(int Imax, int Iils, tData data, boolean verbose) {

        tSolution solut_partial = Solution_init(data);
        tSolution solut_crnt = Solution_init(data);
        tSolution solut_best = Solution_init(data);

        for (int i = 0; i < Imax; ++i) {
            int aux = data.rnd[data.rnd_index++];

            double alpha = R_table(aux);

            if (verbose) {
                System.out.printf("[+] Search %d%n", i + 1);
                System.out.printf("\t[+] Constructing..%n");
            }


            solut_crnt.s = construct(alpha, data);
            update_subseq_info_matrix(solut_crnt, data);

            Solution_cpy(solut_crnt, solut_partial, data);
            if (verbose) {
                System.out.printf("\t[+] Looking for the best Neighbor..%n");
                System.out.printf(Locale.ROOT, "\t    Construction Cost: %.3f%n", solut_partial.cost);
            }

            int iterILS = 0;
            while (iterILS < Iils) {
                RVND(solut_crnt, data);
                if (solut_crnt.cost < solut_partial.cost - DBL_EPSILON) {
                    Solution_cpy(solut_crnt, solut_partial, data);
                    iterILS = 0;
                }

                solut_crnt.s = perturb(solut_partial, data);
                update_subseq_info_matrix(solut_crnt, data);
                iterILS++;
            }

            if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
                Solution_cpy(solut_partial, solut_best, data);
            }

            if (verbose) {
                System.out.println("\tCurrent best cost: " + fmt(solut_best.cost));

                StringBuilder sb = new StringBuilder("SOLUCAO: ");
                for (int k = 0; k < solut_best.s.length; k++) {
                    sb.append(solut_best.s[k]).append(" ");
                }
                System.out.println(sb.toString());
            }

        }
        System.out.printf(Locale.ROOT, "COST: %.2f%n", solut_best.cost);
    }

    // mimics the default formatting of std::cout for doubles (6 significant digits)
    static String fmt(double v) {
        BigDecimal bd = new BigDecimal(v).round(new MathContext(6));
        double d = bd.doubleValue();
        double a = Math.abs(d);
        if (d != 0 && (a < 1e-4 || a >= 1e6)) {
            String s = String.format(Locale.ROOT, "%.5e", d);
            String[] parts = s.split("e");
            String mant = parts[0];
            if (mant.contains(".")) {
                mant = mant.replaceAll("0+$", "").replaceAll("\\.$", "");
            }
            return mant + "e" + parts[1];
        }
        return bd.stripTrailingZeros().toPlainString();
    }

    public static void main(String[] argv) {
        int Imax = 10;
        int Iils;
        boolean verbose = false;

        for (int i = 0; i < argv.length; i++) {
            if (argv[i].equals("-v") || argv[i].equals("--verbose")) {
                verbose = true;
            }
        }

        tData data = new tData();


        List<Integer> rnd = new ArrayList<>();
        double[][][] holder = new double[1][][];

        data.dimen = Data.loadData(holder, rnd);
        data.cost = holder[0];
        data.rnd = new int[rnd.size()];
        for (int i = 0; i < rnd.size(); i++) data.rnd[i] = rnd.get(i);
        data.rnd_index = 0;

        tSolution solut = Solution_init(data);

        Random rand = new Random(System.nanoTime());

        Iils = data.dimen < 100 ? data.dimen : 100;
        long t1 = System.nanoTime();
        ThreadMXBean bean = ManagementFactory.getThreadMXBean();
        long start = bean.getCurrentThreadCpuTime();
        GILS_RVND(Imax, Iils, data, verbose);
        double cpu_time = (double) (bean.getCurrentThreadCpuTime() - start) / 1e9;
        long t2 = System.nanoTime();
        long duration = (t2 - t1) / 1_000_000L;

        double res = (double) duration / 10e2;
        System.out.println("TIME: " + fmt(cpu_time));
        System.out.println("wall clock (s): " + fmt(res));

        System.out.println("Tamanho RND " + rnd.size());

        System.exit(0);
    }
}
