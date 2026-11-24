#define FLAT

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;

namespace MlpSolver
{
    public class TData
    {
        public double[][] cost;
        public int dimen;
        public List<int> rnd;
        public uint rnd_index;
    }

    public struct TSeqInfo
    {
        public double T;
        public double C;
        public double W;
    }

    public enum Neighborhood
    {
        REINSERTION = 1,
        OR_OPT_2,
        OR_OPT_3,
        SWAP,
        TWO_OPT
    }

    public class TSolution
    {
        public List<int> s;
#if FLAT
        public TSeqInfo[] seq;
#elif MATRIX
        public TSeqInfo[][] seq;
#endif
        public double cost;
    }

    public static class Program
    {
        private const int TABLE_SZ = 26;
        private const int DBL_SZ = 8;
        private const int INT_SZ = 4;

        // C's DBL_EPSILON (C#'s double.Epsilon is the smallest denormal, not the same thing)
        private const double DBL_EPSILON = 2.2204460492503131e-16;
        private const double DBL_MAX = double.MaxValue;

        private static int to_1D(int x, int y, int size)
        {
            int a = x * (size + 1) + y;

            return a;
        }

        private static TSolution Solution_init(TData data)
        {
            TSolution solut = new TSolution();
            solut.s = new List<int>(new int[data.dimen + 1]);

#if FLAT
            solut.seq = new TSeqInfo[(data.dimen + 1) * (data.dimen + 1)];
#elif MATRIX
            solut.seq = new TSeqInfo[data.dimen + 1][];
            for (int i = 0; i < data.dimen + 1; i++)
            {
                solut.seq[i] = new TSeqInfo[data.dimen + 1];
            }
#endif
            solut.cost = DBL_MAX;

            return solut;
        }

        private static void Solution_cpy(TSolution src, TSolution tgt, TData data)
        {
            tgt.s = new List<int>(src.s);
            tgt.cost = src.cost;
        }

        private static readonly double[] R_table_values = {
            0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25
        };

        private static double R_table(int i)
        {
            return R_table_values[i];
        }

        private static void print_s(List<int> s)
        {
            for (int i = 0; i < s.Count; i++)
                Console.Write((s[i] + 1) + " ");
            Console.WriteLine();
        }

        private static int partition(List<int> arr, int left, int right, TData data, int r)
        {
            int pivot = arr[right];
            int i = left - 1;
            for (int j = left; j < right; j++)
            {
                if (data.cost[r][arr[j]] < data.cost[r][pivot])
                {
                    i++;
                    swap(arr, i, j);
                }
            }
            swap(arr, i + 1, right);
            return i + 1;
        }

        private static void quicksort(List<int> arr, int left, int right, TData data, int r)
        {
            if (left < right)
            {
                int pivot = partition(arr, left, right, data, r);
                quicksort(arr, left, pivot - 1, data, r);
                quicksort(arr, pivot + 1, right, data, r);
            }
        }

        private static void sort(List<int> arr, int r, TData data)
        {
            quicksort(arr, 0, arr.Count - 1, data, r);
        }

        private static List<int> construct(double alpha, TData data)
        {
            List<int> s = new List<int>(data.dimen + 1);
            s.Add(0);
            List<int> cL = new List<int>(new int[data.dimen - 1]);
            for (int i = 1; i < data.dimen; ++i)
            {
                cL[i - 1] = i;
            }

            int r = 0;
            while (cL.Count != 0)
            {
                sort(cL, r, data);

                int index = data.rnd[(int)data.rnd_index++];

                int c = cL[index];
                s.Add(c);
                r = c;
                cL.RemoveAt(index);
            }

            s.Add(0);

            return s;
        }

        private static void swap(List<int> vec, int i, int j)
        {
            int tmp = vec[i];
            vec[i] = vec[j];
            vec[j] = tmp;
        }

        private static void reverse(List<int> vec, int i, int j)
        {
            vec.Reverse(i, j - i + 1);
        }

        private static void reinsert(List<int> vec, int i, int j, int pos)
        {
            List<int> seq = vec.GetRange(i, j - i + 1);
            if (pos < i)
            {
                vec.RemoveRange(i, j - i + 1);
                vec.InsertRange(pos, seq);
            }
            else
            {
                vec.InsertRange(pos, seq);
                vec.RemoveRange(i, j - i + 1);
            }
        }

#if FLAT

        private static void update_subseq_info_matrix(TSolution solut, TData data, int index = 0)
        {
            int i, j, j_prev, k;
            int from = index;
            bool t;
            for (i = 0; i < data.dimen + 1; i++)
            {
                k = 1 - i - (i == 0 ? 1 : 0);
                t = i == from;

                solut.seq[to_1D(i, i, data.dimen)].T = 0.0;
                solut.seq[to_1D(i, i, data.dimen)].C = 0.0;
                solut.seq[to_1D(i, i, data.dimen)].W = (i == 0 ? 0.0 : 1.0);

                for (j = i + 1; j < data.dimen + 1; j++)
                {
                    j_prev = j - 1;

                    solut.seq[to_1D(i, j, data.dimen)].T = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[to_1D(i, j_prev, data.dimen)].T;
                    solut.seq[to_1D(i, j, data.dimen)].C = solut.seq[to_1D(i, j, data.dimen)].T + solut.seq[to_1D(i, j_prev, data.dimen)].C;
                    solut.seq[to_1D(i, j, data.dimen)].W = j + k;
                }
                from += t ? 1 : 0;
            }

            solut.cost = solut.seq[to_1D(0, data.dimen, data.dimen)].C;
        }

        private static bool search_swap(TSolution solut, TData data)
        {
            double cost_new,
                cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
            double cost_best = DBL_MAX;
            int i, j, j_prev, j_next, i_prev, i_next;
            int I = 0;
            int J = 0;

            for (i = 1; i < data.dimen - 1; ++i)
            {
                i_prev = i - 1;
                i_next = i + 1;

                //consecutive nodes
                cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[i_next]];
                cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, i_next, data.dimen)].T + data.cost[solut.s[i]][solut.s[i_next + 1]];

                cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C +           //       1st subseq
                solut.seq[to_1D(i, i_next, data.dimen)].W * (cost_concat_1) + data.cost[solut.s[i_next]][solut.s[i]] +           // concat 2nd subseq
                solut.seq[to_1D(i_next + 1, data.dimen, data.dimen)].W * (cost_concat_2) + solut.seq[to_1D(i_next + 1, data.dimen, data.dimen)].C;   // concat 3rd subseq

                if (cost_new < cost_best)
                {
                    cost_best = cost_new - DBL_EPSILON;
                    I = i;
                    J = i_next;
                }

                for (j = i_next + 1; j < data.dimen; ++j)
                {
                    j_next = j + 1;
                    j_prev = j - 1;

                    cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[j]];
                    cost_concat_2 = cost_concat_1 + data.cost[solut.s[j]][solut.s[i_next]];
                    cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i_next, j_prev, data.dimen)].T + data.cost[solut.s[j_prev]][solut.s[i]];
                    cost_concat_4 = cost_concat_3 + data.cost[solut.s[i]][solut.s[j_next]];

                    cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C +      // 1st subseq
                    cost_concat_1 +                                                             // concat 2nd subseq (single node)
                    solut.seq[to_1D(i_next, j_prev, data.dimen)].W * cost_concat_2 + solut.seq[to_1D(i_next, j_prev, data.dimen)].C +      // concat 3rd subseq
                    cost_concat_3 +                                                             // concat 4th subseq (single node)
                    solut.seq[to_1D(j_next, data.dimen, data.dimen)].W * cost_concat_4 + solut.seq[to_1D(j_next, data.dimen, data.dimen)].C;   // concat 5th subseq

                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if (cost_best < solut.cost - DBL_EPSILON)
            {
                swap(solut.s, I, J);
                update_subseq_info_matrix(solut, data, I);
                return true;
            }

            return false;
        }

        private static bool search_two_opt(TSolution solut, TData data)
        {
            double cost_new,
                cost_concat_1, cost_concat_2;
            double cost_best = DBL_MAX; // cost_l1, cost_l2;
            double rev_seq_cost;
            int i, j, i_prev, j_next;
            int I = 0;
            int J = 0;

            for (i = 1; i < data.dimen - 1; ++i)
            {
                i_prev = i - 1;

                rev_seq_cost = solut.seq[to_1D(i, i + 1, data.dimen)].T;
                for (j = i + 2; j < data.dimen; ++j)
                {
                    j_next = j + 1;

                    rev_seq_cost += data.cost[solut.s[j - 1]][solut.s[j]] * (solut.seq[to_1D(i, j, data.dimen)].W - 1.0);

                    cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[j]][solut.s[i_prev]];
                    cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, data.dimen)].T + data.cost[solut.s[j_next]][solut.s[i]];

                    cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C +   //  1st subseq
                            solut.seq[to_1D(i, j, data.dimen)].W * cost_concat_1 + rev_seq_cost +   // concat 2nd subseq (reversed seq)
                            solut.seq[to_1D(j_next, data.dimen, data.dimen)].W * cost_concat_2 + solut.seq[to_1D(j_next, data.dimen, data.dimen)].C;      // concat 3rd subseq

                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if (cost_best < solut.cost - DBL_EPSILON)
            {
                reverse(solut.s, I, J);
                update_subseq_info_matrix(solut, data);
                return true;
            }

            return false;
        }

        private static bool search_reinsertion(TSolution solut, TData data, int opt)
        {
            double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
            double cost_best = DBL_MAX; //, cost_l1, cost_l2, cost_l3;
            int i, j, k, k_next, i_prev, j_next;
            int I = 0;
            int J = 0;
            int POS = 0;

            for (i = 1, j = opt + i - 1; i < data.dimen - opt + 1; ++i, ++j)
            {
                j_next = j + 1;
                i_prev = i - 1;

                //k -> edges
                for (k = 0; k < i_prev; ++k)
                {
                    k_next = k + 1;

                    cost_concat_1 = solut.seq[to_1D(0, k, data.dimen)].T + data.cost[solut.s[k]][solut.s[i]];
                    cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, data.dimen)].T + data.cost[solut.s[j]][solut.s[k_next]];
                    cost_concat_3 = cost_concat_2 + solut.seq[to_1D(k_next, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]];

                    cost_new = solut.seq[to_1D(0, k, data.dimen)].C +   //       1st subseq
                        solut.seq[to_1D(i, j, data.dimen)].W * cost_concat_1 + solut.seq[to_1D(i, j, data.dimen)].C +   //  concat 2nd subseq (reinserted seq)
                        solut.seq[to_1D(k_next, i_prev, data.dimen)].W * cost_concat_2 + solut.seq[to_1D(k_next, i_prev, data.dimen)].C +   //  concat 3rd subseq
                        solut.seq[to_1D(j_next, data.dimen, data.dimen)].W * cost_concat_3 + solut.seq[to_1D(j_next, data.dimen, data.dimen)].C;       // concat 4th subseq

                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }

                for (k = i + opt; k < data.dimen; ++k)
                {
                    k_next = k + 1;

                    cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]];
                    cost_concat_2 = cost_concat_1 + solut.seq[to_1D(j_next, k, data.dimen)].T + data.cost[solut.s[k]][solut.s[i]];
                    cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i, j, data.dimen)].T + data.cost[solut.s[j]][solut.s[k_next]];

                    cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C +   //       1st subseq
                        solut.seq[to_1D(j_next, k, data.dimen)].W * cost_concat_1 + solut.seq[to_1D(j_next, k, data.dimen)].C +   // concat 2nd subseq
                        solut.seq[to_1D(i, j, data.dimen)].W * cost_concat_2 + solut.seq[to_1D(i, j, data.dimen)].C +   // concat 3rd subseq (reinserted seq)
                        solut.seq[to_1D(k_next, data.dimen, data.dimen)].W * cost_concat_3 + solut.seq[to_1D(k_next, data.dimen, data.dimen)].C;       // concat 4th subseq


                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }
            }

            if (cost_best < solut.cost - DBL_EPSILON)
            {
                reinsert(solut.s, I, J, POS + 1);
                update_subseq_info_matrix(solut, data, I < POS + 1 ? I : POS + 1);
                return true;
            }

            return false;
        }

#elif MATRIX

        private static void update_subseq_info_matrix(TSolution solut, TData data, int index = 0)
        {
            int i, j, j_prev, k;
            int from = index;
            bool t;
            for (i = 0; i < data.dimen + 1; i++)
            {
                k = 1 - i - (i == 0 ? 1 : 0);
                t = i == from;

                solut.seq[i][i].T = 0.0;
                solut.seq[i][i].C = 0.0;
                solut.seq[i][i].W = (i == 0 ? 0.0 : 1.0);

                for (j = i + 1; j < data.dimen + 1; j++)
                {
                    j_prev = j - 1;

                    solut.seq[i][j].T = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev].T;
                    solut.seq[i][j].C = solut.seq[i][j].T + solut.seq[i][j_prev].C;
                    solut.seq[i][j].W = j + k;

                }
                from += t ? 1 : 0;
            }

            solut.cost = solut.seq[0][data.dimen].C;
        }

        private static bool search_swap(TSolution solut, TData data)
        {
            double cost_new,
                cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
            double cost_best = DBL_MAX;
            int i, j, j_prev, j_next, i_prev, i_next;
            int I = 0;
            int J = 0;

            for (i = 1; i < data.dimen - 1; ++i)
            {
                i_prev = i - 1;
                i_next = i + 1;

                //consecutive nodes
                cost_concat_1 = solut.seq[0][i_prev].T + data.cost[solut.s[i_prev]][solut.s[i_next]];
                cost_concat_2 = cost_concat_1 + solut.seq[i][i_next].T + data.cost[solut.s[i]][solut.s[i_next + 1]];

                cost_new = solut.seq[0][i_prev].C +           //       1st subseq
                solut.seq[i][i_next].W * (cost_concat_1) + data.cost[solut.s[i_next]][solut.s[i]] +           // concat 2nd subseq
                solut.seq[i_next + 1][data.dimen].W * (cost_concat_2) + solut.seq[i_next + 1][data.dimen].C;   // concat 3rd subseq

                if (cost_new < cost_best)
                {
                    cost_best = cost_new - DBL_EPSILON;
                    I = i;
                    J = i_next;
                }

                for (j = i_next + 1; j < data.dimen; ++j)
                {
                    j_next = j + 1;
                    j_prev = j - 1;

                    cost_concat_1 = solut.seq[0][i_prev].T + data.cost[solut.s[i_prev]][solut.s[j]];
                    cost_concat_2 = cost_concat_1 + data.cost[solut.s[j]][solut.s[i_next]];
                    cost_concat_3 = cost_concat_2 + solut.seq[i_next][j_prev].T + data.cost[solut.s[j_prev]][solut.s[i]];
                    cost_concat_4 = cost_concat_3 + data.cost[solut.s[i]][solut.s[j_next]];

                    cost_new = solut.seq[0][i_prev].C +      // 1st subseq
                    cost_concat_1 +                                                             // concat 2nd subseq (single node)
                    solut.seq[i_next][j_prev].W * cost_concat_2 + solut.seq[i_next][j_prev].C +      // concat 3rd subseq
                    cost_concat_3 +                                                             // concat 4th subseq (single node)
                    solut.seq[j_next][data.dimen].W * cost_concat_4 + solut.seq[j_next][data.dimen].C;   // concat 5th subseq

                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if (cost_best < solut.cost - DBL_EPSILON)
            {
                swap(solut.s, I, J);
                update_subseq_info_matrix(solut, data, I);
                return true;
            }

            return false;
        }

        private static bool search_two_opt(TSolution solut, TData data)
        {
            double cost_new,
                cost_concat_1, cost_concat_2;
            double cost_best = DBL_MAX; // cost_l1, cost_l2;
            double rev_seq_cost;
            int i, j, i_prev, j_next;
            int I = 0;
            int J = 0;

            for (i = 1; i < data.dimen - 1; ++i)
            {
                i_prev = i - 1;

                rev_seq_cost = solut.seq[i][i + 1].T;
                for (j = i + 2; j < data.dimen; ++j)
                {
                    j_next = j + 1;


                    rev_seq_cost += data.cost[solut.s[j - 1]][solut.s[j]] * (solut.seq[i][j].W - 1.0);

                    cost_concat_1 = solut.seq[0][i_prev].T + data.cost[solut.s[j]][solut.s[i_prev]];
                    cost_concat_2 = cost_concat_1 + solut.seq[i][j].T + data.cost[solut.s[j_next]][solut.s[i]];

                    cost_new = solut.seq[0][i_prev].C +   //  1st subseq
                        solut.seq[i][j].W * cost_concat_1 + rev_seq_cost +   // concat 2nd subseq (reversed seq)
                        solut.seq[j_next][data.dimen].W * cost_concat_2 + solut.seq[j_next][data.dimen].C;      // concat 3rd subseq


                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                    }
                }
            }

            if (cost_best < solut.cost - DBL_EPSILON)
            {
                reverse(solut.s, I, J);
                update_subseq_info_matrix(solut, data);
                return true;
            }

            return false;
        }

        private static bool search_reinsertion(TSolution solut, TData data, int opt)
        {
            double cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
            double cost_best = DBL_MAX; //, cost_l1, cost_l2, cost_l3;
            int i, j, k, k_next, i_prev, j_next;
            int I = 0;
            int J = 0;
            int POS = 0;

            for (i = 1, j = opt + i - 1; i < data.dimen - opt + 1; ++i, ++j)
            {
                j_next = j + 1;
                i_prev = i - 1;

                //k -> edges
                for (k = 0; k < i_prev; ++k)
                {
                    k_next = k + 1;

                    cost_concat_1 = solut.seq[0][k].T + data.cost[solut.s[k]][solut.s[i]];
                    cost_concat_2 = cost_concat_1 + solut.seq[i][j].T + data.cost[solut.s[j]][solut.s[k_next]];
                    cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev].T + data.cost[solut.s[i_prev]][solut.s[j_next]];

                    cost_new = solut.seq[0][k].C +   //       1st subseq
                        solut.seq[i][j].W * cost_concat_1 + solut.seq[i][j].C +   //  concat 2nd subseq (reinserted seq)
                        solut.seq[k_next][i_prev].W * cost_concat_2 + solut.seq[k_next][i_prev].C +   //  concat 3rd subseq
                        solut.seq[j_next][data.dimen].W * cost_concat_3 + solut.seq[j_next][data.dimen].C;       // concat 4th subseq

                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }

                for (k = i + opt; k < data.dimen; ++k)
                {
                    k_next = k + 1;

                    cost_concat_1 = solut.seq[0][i_prev].T + data.cost[solut.s[i_prev]][solut.s[j_next]];
                    cost_concat_2 = cost_concat_1 + solut.seq[j_next][k].T + data.cost[solut.s[k]][solut.s[i]];
                    cost_concat_3 = cost_concat_2 + solut.seq[i][j].T + data.cost[solut.s[j]][solut.s[k_next]];

                    cost_new = solut.seq[0][i_prev].C +   //       1st subseq
                            solut.seq[j_next][k].W * cost_concat_1 + solut.seq[j_next][k].C +   // concat 2nd subseq
                            solut.seq[i][j].W * cost_concat_2 + solut.seq[i][j].C +   // concat 3rd subseq (reinserted seq)
                            solut.seq[k_next][data.dimen].W * cost_concat_3 + solut.seq[k_next][data.dimen].C;       // concat 4th subseq

                    if (cost_new < cost_best)
                    {
                        cost_best = cost_new - DBL_EPSILON;
                        I = i;
                        J = j;
                        POS = k;
                    }
                }
            }

            if (cost_best < solut.cost - DBL_EPSILON)
            {
                reinsert(solut.s, I, J, POS + 1);
                update_subseq_info_matrix(solut, data, I < POS + 1 ? I : POS + 1);
                return true;
            }

            return false;
        }

#endif

        private static void RVND(TSolution solut, TData data)
        {
            List<Neighborhood> neighbd_list = new List<Neighborhood> { Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION, Neighborhood.OR_OPT_2, Neighborhood.OR_OPT_3 };
            int index;
            Neighborhood neighbd;
            bool improve_flag;

            while (neighbd_list.Count != 0)
            {
                index = data.rnd[(int)data.rnd_index++];
                neighbd = neighbd_list[index];
                improve_flag = false;

                switch (neighbd)
                {
                    case Neighborhood.REINSERTION:
                        improve_flag = search_reinsertion(solut, data, (int)Neighborhood.REINSERTION);
                        break;
                    case Neighborhood.OR_OPT_2:
                        improve_flag = search_reinsertion(solut, data, (int)Neighborhood.OR_OPT_2);
                        break;
                    case Neighborhood.OR_OPT_3:
                        improve_flag = search_reinsertion(solut, data, (int)Neighborhood.OR_OPT_3);
                        break;
                    case Neighborhood.SWAP:
                        improve_flag = search_swap(solut, data);
                        break;
                    case Neighborhood.TWO_OPT:
                        improve_flag = search_two_opt(solut, data);
                        break;
                }
                if (improve_flag)
                {
                    neighbd_list = new List<Neighborhood> { Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION, Neighborhood.OR_OPT_2, Neighborhood.OR_OPT_3 };
                }
                else
                {
                    neighbd_list.RemoveAt(index);
                }
            }
        }

        private static List<int> perturb(TSolution solut, TData data)
        {
            List<int> s = new List<int>(solut.s);
            int A_start = 1;
            int A_end = 1;
            int B_start = 1;
            int B_end = 1;

            int size_max = (data.dimen + 1) / 10;
            size_max = size_max >= 2 ? size_max : 2;
            int size_min = 2;
            while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end))
            {
                A_start = data.rnd[(int)data.rnd_index++];
                A_end = A_start + data.rnd[(int)data.rnd_index++];

                B_start = data.rnd[(int)data.rnd_index++];
                B_end = B_start + data.rnd[(int)data.rnd_index++];
            }

            if (A_start < B_start)
            {
                reinsert(s, B_start, B_end - 1, A_end);
                reinsert(s, A_start, A_end - 1, B_end);
            }
            else
            {
                reinsert(s, A_start, A_end - 1, B_end);
                reinsert(s, B_start, B_end - 1, A_end);
            }

            return s;
        }

        private static void GILS_RVND(int Imax, int Iils, TData data, bool verbose)
        {
            TSolution solut_partial = Solution_init(data);
            TSolution solut_crnt = Solution_init(data);
            TSolution solut_best = Solution_init(data);

            for (int i = 0; i < Imax; ++i)
            {
                int aux = data.rnd[(int)data.rnd_index++];

                double alpha = R_table(aux);

                if (verbose)
                {
                    Console.Write(string.Format(CultureInfo.InvariantCulture, "[+] Search {0}\n", i + 1));
                    Console.Write("\t[+] Constructing..\n");
                }

                solut_crnt.s = construct(alpha, data);
                update_subseq_info_matrix(solut_crnt, data);

                Solution_cpy(solut_crnt, solut_partial, data);
                if (verbose)
                {
                    Console.Write("\t[+] Looking for the best Neighbor..\n");
                    Console.Write(string.Format(CultureInfo.InvariantCulture, "\t    Construction Cost: {0:F3}\n", solut_partial.cost));
                }

                int iterILS = 0;
                while (iterILS < Iils)
                {
                    RVND(solut_crnt, data);
                    if (solut_crnt.cost < solut_partial.cost - DBL_EPSILON)
                    {
                        Solution_cpy(solut_crnt, solut_partial, data);
                        iterILS = 0;
                    }

                    solut_crnt.s = perturb(solut_partial, data);
                    update_subseq_info_matrix(solut_crnt, data);
                    iterILS++;
                }

                if (solut_partial.cost < solut_best.cost - DBL_EPSILON)
                {
                    Solution_cpy(solut_partial, solut_best, data);
                }

                if (verbose)
                {
                    Console.WriteLine("\tCurrent best cost: " + Ostream(solut_best.cost));

                    Console.Write("SOLUCAO: ");
                    for (int k = 0; k < solut_best.s.Count; k++)
                    {
                        Console.Write(solut_best.s[k] + " ");
                    }
                    Console.WriteLine();
                }
            }
            Console.Write(string.Format(CultureInfo.InvariantCulture, "COST: {0:F2}\n", solut_best.cost));
        }

        // Mimics the default formatting of std::ostream for doubles (6 significant digits)
        private static string Ostream(double value)
        {
            return value.ToString("G6", CultureInfo.InvariantCulture);
        }

        public static int Main(string[] argv)
        {
            int Imax = 10;
            int Iils;
            bool verbose = false;

            for (int i = 0; i < argv.Length; i++)
            {
                if (argv[i] == "-v" || argv[i] == "--verbose")
                {
                    verbose = true;
                }
            }

            TData data = new TData();

            List<int> rnd = new List<int>();

            data.dimen = Data.LoadData(out data.cost, rnd);
            data.rnd = rnd;
            data.rnd_index = 0;

            TSolution solut = Solution_init(data);

            Iils = data.dimen < 100 ? data.dimen : 100;
            Stopwatch t1 = Stopwatch.StartNew();
            TimeSpan start = Process.GetCurrentProcess().TotalProcessorTime;
            GILS_RVND(Imax, Iils, data, verbose);
            double cpu_time = (Process.GetCurrentProcess().TotalProcessorTime - start).TotalSeconds;
            t1.Stop();
            long duration = (long)t1.Elapsed.TotalMilliseconds;

            double res = (double)duration / 10e2;
            Console.WriteLine("TIME: " + Ostream(cpu_time));
            Console.WriteLine("wall clock (s): " + Ostream(res));

            Console.WriteLine("Tamanho RND " + rnd.Count);

            Environment.Exit(0);

            return 0;
        }
    }
}
