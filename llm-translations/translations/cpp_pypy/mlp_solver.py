#!/usr/bin/env python3
"""GILS-RVND solver for the Minimum Latency Problem (translated from main.cpp)."""

import sys
import time
from enum import IntEnum

sys.setrecursionlimit(100000)

EPS = sys.float_info.epsilon
DBL_MAX = sys.float_info.max


class Neighborhood(IntEnum):
    REINSERTION = 1
    OR_OPT_2 = 2
    OR_OPT_3 = 3
    SWAP = 4
    TWO_OPT = 5


class Data:
    __slots__ = ("cost", "dimen", "rnd", "rnd_index")


class Solution:
    __slots__ = ("s", "seqT", "seqC", "seqW", "cost")


def R_table(i):
    table = (0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
              0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25)
    return table[i]


def print_s(s):
    print(" ".join(str(x + 1) for x in s))


def to_1D(x, y, size):
    return x * (size + 1) + y


def matrix_allocate(sz):
    return [[0.0] * sz for _ in range(sz)]


def loadData():
    rnd = []
    with open("../distance_matrix", "r") as f:
        first_line = f.readline()
        dimension = int(first_line.strip())
        matrix = matrix_allocate(dimension)

        flag = False
        rnd_size = -1
        i = 0
        for raw_line in f:
            line = raw_line.rstrip("\n").rstrip("\r")
            j = i + 1
            while " " in line:
                pos = line.find(" ")
                value = float(line[:pos])
                line = line[pos + 1:]
                matrix[i][j] = value
                matrix[j][i] = value
                j += 1

            i += 1

            if flag:
                if rnd_size < 0:
                    rnd_size = int(line)
                else:
                    rnd.append(int(line))

            if line == "RND":
                flag = True

    return dimension, matrix, rnd


def Solution_init(data):
    solut = Solution()
    solut.s = [0] * (data.dimen + 1)
    n = (data.dimen + 1) * (data.dimen + 1)
    solut.seqT = [0.0] * n
    solut.seqC = [0.0] * n
    solut.seqW = [0.0] * n
    solut.cost = DBL_MAX
    return solut


def Solution_cpy(src, tgt):
    tgt.s = list(src.s)
    tgt.cost = src.cost


def partition(arr, left, right, cost, r):
    pivot = arr[right]
    row = cost[r]
    i = left - 1
    for j in range(left, right):
        if row[arr[j]] < row[pivot]:
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1


def quicksort(arr, left, right, cost, r):
    if left < right:
        p = partition(arr, left, right, cost, r)
        quicksort(arr, left, p - 1, cost, r)
        quicksort(arr, p + 1, right, cost, r)


def sort(arr, r, data):
    quicksort(arr, 0, len(arr) - 1, data.cost, r)


def construct(alpha, data):
    s = [0]
    cL = list(range(1, data.dimen))

    r = 0
    while cL:
        sort(cL, r, data)

        index = data.rnd[data.rnd_index]
        data.rnd_index += 1

        c = cL[index]
        s.append(c)
        r = c
        del cL[index]

    s.append(0)
    return s


def reinsert(vec, i, j, pos):
    seq = vec[i:j + 1]
    if pos < i:
        del vec[i:j + 1]
        vec[pos:pos] = seq
    else:
        vec[pos:pos] = seq
        del vec[i:j + 1]


def update_subseq_info_matrix(solut, data, index=0):
    dimen = data.dimen
    stride = dimen + 1
    cost = data.cost
    s = solut.s
    seqT = solut.seqT
    seqC = solut.seqC
    seqW = solut.seqW

    from_ = index
    for i in range(dimen + 1):
        k = 1 - i - (1 if i == 0 else 0)
        t = (i == from_)

        idx_ii = i * stride + i
        seqT[idx_ii] = 0.0
        seqC[idx_ii] = 0.0
        seqW[idx_ii] = 0.0 if i == 0 else 1.0

        base = i * stride
        for j in range(i + 1, dimen + 1):
            j_prev = j - 1
            idx_ij = base + j
            idx_ijprev = base + j_prev

            seqT[idx_ij] = cost[s[j_prev]][s[j]] + seqT[idx_ijprev]
            seqC[idx_ij] = seqT[idx_ij] + seqC[idx_ijprev]
            seqW[idx_ij] = j + k

        if t:
            from_ += 1

    solut.cost = seqC[dimen]


def search_swap(solut, data):
    dimen = data.dimen
    stride = dimen + 1
    cost = data.cost
    s = solut.s
    seqT = solut.seqT
    seqC = solut.seqC
    seqW = solut.seqW

    cost_best = DBL_MAX
    I = J = -1

    for i in range(1, dimen - 1):
        i_prev = i - 1
        i_next = i + 1

        idx_0_iprev = i_prev
        idx_i_inext = i * stride + i_next
        idx_inextp1_dimen = (i_next + 1) * stride + dimen

        cost_concat_1 = seqT[idx_0_iprev] + cost[s[i_prev]][s[i_next]]
        cost_concat_2 = cost_concat_1 + seqT[idx_i_inext] + cost[s[i]][s[i_next + 1]]

        cost_new = (seqC[idx_0_iprev]
                    + seqW[idx_i_inext] * cost_concat_1 + cost[s[i_next]][s[i]]
                    + seqW[idx_inextp1_dimen] * cost_concat_2 + seqC[idx_inextp1_dimen])

        if cost_new < cost_best:
            cost_best = cost_new - EPS
            I = i
            J = i_next

        for j in range(i_next + 1, dimen):
            j_next = j + 1
            j_prev = j - 1

            idx_inext_jprev = i_next * stride + j_prev
            idx_jnext_dimen = j_next * stride + dimen

            cost_concat_1 = seqT[idx_0_iprev] + cost[s[i_prev]][s[j]]
            cost_concat_2 = cost_concat_1 + cost[s[j]][s[i_next]]
            cost_concat_3 = cost_concat_2 + seqT[idx_inext_jprev] + cost[s[j_prev]][s[i]]
            cost_concat_4 = cost_concat_3 + cost[s[i]][s[j_next]]

            cost_new = (seqC[idx_0_iprev]
                        + cost_concat_1
                        + seqW[idx_inext_jprev] * cost_concat_2 + seqC[idx_inext_jprev]
                        + cost_concat_3
                        + seqW[idx_jnext_dimen] * cost_concat_4 + seqC[idx_jnext_dimen])

            if cost_new < cost_best:
                cost_best = cost_new - EPS
                I = i
                J = j

    if cost_best < solut.cost - EPS:
        s[I], s[J] = s[J], s[I]
        update_subseq_info_matrix(solut, data, I)
        return True

    return False


def search_two_opt(solut, data):
    dimen = data.dimen
    stride = dimen + 1
    cost = data.cost
    s = solut.s
    seqT = solut.seqT
    seqC = solut.seqC
    seqW = solut.seqW

    cost_best = DBL_MAX
    I = J = -1

    for i in range(1, dimen - 1):
        i_prev = i - 1
        idx_0_iprev = i_prev

        rev_seq_cost = seqT[i * stride + i + 1]
        for j in range(i + 2, dimen):
            j_next = j + 1
            idx_i_j = i * stride + j

            rev_seq_cost += cost[s[j - 1]][s[j]] * (seqW[idx_i_j] - 1.0)

            idx_jnext_dimen = j_next * stride + dimen

            cost_concat_1 = seqT[idx_0_iprev] + cost[s[j]][s[i_prev]]
            cost_concat_2 = cost_concat_1 + seqT[idx_i_j] + cost[s[j_next]][s[i]]

            cost_new = (seqC[idx_0_iprev]
                        + seqW[idx_i_j] * cost_concat_1 + rev_seq_cost
                        + seqW[idx_jnext_dimen] * cost_concat_2 + seqC[idx_jnext_dimen])

            if cost_new < cost_best:
                cost_best = cost_new - EPS
                I = i
                J = j

    if cost_best < solut.cost - EPS:
        s[I:J + 1] = s[I:J + 1][::-1]
        update_subseq_info_matrix(solut, data)
        return True

    return False


def search_reinsertion(solut, data, opt):
    dimen = data.dimen
    stride = dimen + 1
    cost = data.cost
    s = solut.s
    seqT = solut.seqT
    seqC = solut.seqC
    seqW = solut.seqW

    cost_best = DBL_MAX
    I = J = POS = -1

    i = 1
    j = opt + i - 1
    while i < dimen - opt + 1:
        j_next = j + 1
        i_prev = i - 1

        idx_i_j = i * stride + j

        for k in range(0, i_prev):
            k_next = k + 1

            idx_0_k = k
            idx_knext_iprev = k_next * stride + i_prev
            idx_jnext_dimen = j_next * stride + dimen

            cost_concat_1 = seqT[idx_0_k] + cost[s[k]][s[i]]
            cost_concat_2 = cost_concat_1 + seqT[idx_i_j] + cost[s[j]][s[k_next]]
            cost_concat_3 = cost_concat_2 + seqT[idx_knext_iprev] + cost[s[i_prev]][s[j_next]]

            cost_new = (seqC[idx_0_k]
                        + seqW[idx_i_j] * cost_concat_1 + seqC[idx_i_j]
                        + seqW[idx_knext_iprev] * cost_concat_2 + seqC[idx_knext_iprev]
                        + seqW[idx_jnext_dimen] * cost_concat_3 + seqC[idx_jnext_dimen])

            if cost_new < cost_best:
                cost_best = cost_new - EPS
                I = i
                J = j
                POS = k

        for k in range(i + opt, dimen):
            k_next = k + 1

            idx_0_iprev = i_prev
            idx_jnext_k = j_next * stride + k
            idx_knext_dimen = k_next * stride + dimen

            cost_concat_1 = seqT[idx_0_iprev] + cost[s[i_prev]][s[j_next]]
            cost_concat_2 = cost_concat_1 + seqT[idx_jnext_k] + cost[s[k]][s[i]]
            cost_concat_3 = cost_concat_2 + seqT[idx_i_j] + cost[s[j]][s[k_next]]

            cost_new = (seqC[idx_0_iprev]
                        + seqW[idx_jnext_k] * cost_concat_1 + seqC[idx_jnext_k]
                        + seqW[idx_i_j] * cost_concat_2 + seqC[idx_i_j]
                        + seqW[idx_knext_dimen] * cost_concat_3 + seqC[idx_knext_dimen])

            if cost_new < cost_best:
                cost_best = cost_new - EPS
                I = i
                J = j
                POS = k

        i += 1
        j += 1

    if cost_best < solut.cost - EPS:
        reinsert(solut.s, I, J, POS + 1)
        update_subseq_info_matrix(solut, data, I if I < POS + 1 else POS + 1)
        return True

    return False


def RVND(solut, data):
    neighbd_list = [Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION,
                    Neighborhood.OR_OPT_2, Neighborhood.OR_OPT_3]

    while neighbd_list:
        index = data.rnd[data.rnd_index]
        data.rnd_index += 1
        neighbd = neighbd_list[index]
        improve_flag = False

        if neighbd == Neighborhood.REINSERTION:
            improve_flag = search_reinsertion(solut, data, int(Neighborhood.REINSERTION))
        elif neighbd == Neighborhood.OR_OPT_2:
            improve_flag = search_reinsertion(solut, data, int(Neighborhood.OR_OPT_2))
        elif neighbd == Neighborhood.OR_OPT_3:
            improve_flag = search_reinsertion(solut, data, int(Neighborhood.OR_OPT_3))
        elif neighbd == Neighborhood.SWAP:
            improve_flag = search_swap(solut, data)
        elif neighbd == Neighborhood.TWO_OPT:
            improve_flag = search_two_opt(solut, data)

        if improve_flag:
            neighbd_list = [Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION,
                            Neighborhood.OR_OPT_2, Neighborhood.OR_OPT_3]
        else:
            del neighbd_list[index]


def perturb(solut, data):
    s = list(solut.s)
    A_start = 1
    A_end = 1
    B_start = 1
    B_end = 1

    size_max = (data.dimen + 1) // 10
    size_max = size_max if size_max >= 2 else 2
    size_min = 2

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end):
        A_start = data.rnd[data.rnd_index]
        data.rnd_index += 1
        A_end = A_start + data.rnd[data.rnd_index]
        data.rnd_index += 1

        B_start = data.rnd[data.rnd_index]
        data.rnd_index += 1
        B_end = B_start + data.rnd[data.rnd_index]
        data.rnd_index += 1

    if A_start < B_start:
        reinsert(s, B_start, B_end - 1, A_end)
        reinsert(s, A_start, A_end - 1, B_end)
    else:
        reinsert(s, A_start, A_end - 1, B_end)
        reinsert(s, B_start, B_end - 1, A_end)

    return s


def GILS_RVND(Imax, Iils, data, verbose):
    solut_partial = Solution_init(data)
    solut_crnt = Solution_init(data)
    solut_best = Solution_init(data)

    for i in range(Imax):
        aux = data.rnd[data.rnd_index]
        data.rnd_index += 1

        alpha = R_table(aux)

        if verbose:
            print("[+] Search %d" % (i + 1))
            print("\t[+] Constructing..")

        solut_crnt.s = construct(alpha, data)
        update_subseq_info_matrix(solut_crnt, data)

        Solution_cpy(solut_crnt, solut_partial)
        if verbose:
            print("\t[+] Looking for the best Neighbor..")
            print("\t    Construction Cost: %.3lf" % solut_partial.cost)

        iterILS = 0
        while iterILS < Iils:
            RVND(solut_crnt, data)
            if solut_crnt.cost < solut_partial.cost - EPS:
                Solution_cpy(solut_crnt, solut_partial)
                iterILS = 0

            solut_crnt.s = perturb(solut_partial, data)
            update_subseq_info_matrix(solut_crnt, data)
            iterILS += 1

        if solut_partial.cost < solut_best.cost - EPS:
            Solution_cpy(solut_partial, solut_best)

        if verbose:
            print("\tCurrent best cost: %s" % solut_best.cost)

            print("SOLUCAO: " + " ".join(str(x) for x in solut_best.s))

    print("COST: %.2lf" % solut_best.cost)


def main():
    Imax = 10
    verbose = False

    for arg in sys.argv[1:]:
        if arg == "-v" or arg == "--verbose":
            verbose = True

    data = Data()

    rnd = []
    dimen, cost, rnd = loadData()
    data.cost = cost
    data.dimen = dimen
    data.rnd = rnd
    data.rnd_index = 0

    Solution_init(data)

    Iils = data.dimen if data.dimen < 100 else 100

    t1 = time.perf_counter()
    start = time.process_time()
    GILS_RVND(Imax, Iils, data, verbose)
    cpu_time = time.process_time() - start
    t2 = time.perf_counter()
    duration = (t2 - t1) * 1000

    res = duration / 10e2
    print("TIME: %s" % cpu_time)
    print("wall clock (s): %s" % res)

    print("Tamanho RND %d" % len(rnd))

    sys.exit(0)


if __name__ == "__main__":
    main()
