import sys
import time
import random
from enum import IntEnum

from data import load_data

sys.setrecursionlimit(10000)

DBL_EPSILON = sys.float_info.epsilon

R_TABLE = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
           0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]


class Neighborhood(IntEnum):
    REINSERTION = 1
    OR_OPT_2 = 2
    OR_OPT_3 = 3
    SWAP = 4
    TWO_OPT = 5


class TData:
    def __init__(self, cost, dimen, rnd):
        self.cost = cost
        self.dimen = dimen
        self.rnd = rnd
        self.rnd_index = 0


class SeqInfo:
    __slots__ = ("T", "C", "W")

    def __init__(self):
        self.T = 0.0
        self.C = 0.0
        self.W = 0.0


class Solution:
    def __init__(self, data):
        self.s = [0] * (data.dimen + 1)
        self.seq = [SeqInfo() for _ in range((data.dimen + 1) * (data.dimen + 1))]
        self.cost = float("inf")


def to_1D(x, y, size):
    return x * (size + 1) + y


def solution_copy(src, tgt):
    tgt.s = list(src.s)
    tgt.cost = src.cost


def r_table(i):
    return R_TABLE[i]


def print_s(s):
    print(" ".join(str(v + 1) for v in s))


def partition(arr, left, right, data, r):
    pivot = arr[right]
    i = left - 1
    for j in range(left, right):
        if data.cost[r][arr[j]] < data.cost[r][pivot]:
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1


def quicksort(arr, left, right, data, r):
    if left < right:
        pivot = partition(arr, left, right, data, r)
        quicksort(arr, left, pivot - 1, data, r)
        quicksort(arr, pivot + 1, right, data, r)


def sort_cl(arr, r, data):
    quicksort(arr, 0, len(arr) - 1, data, r)


def construct(alpha, data):
    s = [0]
    cL = list(range(1, data.dimen))

    r = 0
    while cL:
        sort_cl(cL, r, data)

        index = data.rnd[data.rnd_index]
        data.rnd_index += 1

        c = cL[index]
        s.append(c)
        r = c
        del cL[index]

    s.append(0)

    return s


def swap(vec, i, j):
    vec[i], vec[j] = vec[j], vec[i]


def reverse(vec, i, j):
    vec[i:j + 1] = vec[i:j + 1][::-1]


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
    from_ = index

    for i in range(dimen + 1):
        k = 1 - i - (1 if i == 0 else 0)
        t = (i == from_)

        ii = to_1D(i, i, dimen)
        solut.seq[ii].T = 0.0
        solut.seq[ii].C = 0.0
        solut.seq[ii].W = 0.0 if i == 0 else 1.0

        for j in range(i + 1, dimen + 1):
            j_prev = j - 1
            idx = to_1D(i, j, dimen)
            idx_prev = to_1D(i, j_prev, dimen)

            solut.seq[idx].T = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[idx_prev].T
            solut.seq[idx].C = solut.seq[idx].T + solut.seq[idx_prev].C
            solut.seq[idx].W = j + k

        from_ += 1 if t else 0

    solut.cost = solut.seq[to_1D(0, dimen, dimen)].C


def search_swap(solut, data):
    dimen = data.dimen
    cost_best = float("inf")
    I = J = None

    for i in range(1, dimen - 1):
        i_prev = i - 1
        i_next = i + 1

        cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, i_next, dimen)].T + data.cost[solut.s[i]][solut.s[i_next + 1]]

        cost_new = (solut.seq[to_1D(0, i_prev, dimen)].C
                    + solut.seq[to_1D(i, i_next, dimen)].W * cost_concat_1 + data.cost[solut.s[i_next]][solut.s[i]]
                    + solut.seq[to_1D(i_next + 1, dimen, dimen)].W * cost_concat_2
                    + solut.seq[to_1D(i_next + 1, dimen, dimen)].C)

        if cost_new < cost_best:
            cost_best = cost_new - DBL_EPSILON
            I, J = i, i_next

        for j in range(i_next + 1, dimen):
            j_next = j + 1
            j_prev = j - 1

            cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + data.cost[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i_next, j_prev, dimen)].T + data.cost[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3 + data.cost[solut.s[i]][solut.s[j_next]]

            cost_new = (solut.seq[to_1D(0, i_prev, dimen)].C
                        + cost_concat_1
                        + solut.seq[to_1D(i_next, j_prev, dimen)].W * cost_concat_2
                        + solut.seq[to_1D(i_next, j_prev, dimen)].C
                        + cost_concat_3
                        + solut.seq[to_1D(j_next, dimen, dimen)].W * cost_concat_4
                        + solut.seq[to_1D(j_next, dimen, dimen)].C)

            if cost_new < cost_best:
                cost_best = cost_new - DBL_EPSILON
                I, J = i, j

    if cost_best < solut.cost - DBL_EPSILON:
        swap(solut.s, I, J)
        update_subseq_info_matrix(solut, data, I)
        return True

    return False


def search_two_opt(solut, data):
    dimen = data.dimen
    cost_best = float("inf")
    I = J = None

    for i in range(1, dimen - 1):
        i_prev = i - 1

        rev_seq_cost = solut.seq[to_1D(i, i + 1, dimen)].T
        for j in range(i + 2, dimen):
            j_next = j + 1

            rev_seq_cost += data.cost[solut.s[j - 1]][solut.s[j]] * (solut.seq[to_1D(i, j, dimen)].W - 1.0)

            cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, dimen)].T + data.cost[solut.s[j_next]][solut.s[i]]

            cost_new = (solut.seq[to_1D(0, i_prev, dimen)].C
                        + solut.seq[to_1D(i, j, dimen)].W * cost_concat_1 + rev_seq_cost
                        + solut.seq[to_1D(j_next, dimen, dimen)].W * cost_concat_2
                        + solut.seq[to_1D(j_next, dimen, dimen)].C)

            if cost_new < cost_best:
                cost_best = cost_new - DBL_EPSILON
                I, J = i, j

    if cost_best < solut.cost - DBL_EPSILON:
        reverse(solut.s, I, J)
        update_subseq_info_matrix(solut, data)
        return True

    return False


def search_reinsertion(solut, data, opt):
    dimen = data.dimen
    cost_best = float("inf")
    I = J = POS = None

    i = 1
    j = opt + i - 1
    while i < dimen - opt + 1:
        j_next = j + 1
        i_prev = i - 1

        for k in range(0, i_prev):
            k_next = k + 1

            cost_concat_1 = solut.seq[to_1D(0, k, dimen)].T + data.cost[solut.s[k]][solut.s[i]]
            cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, dimen)].T + data.cost[solut.s[j]][solut.s[k_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[to_1D(k_next, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]]

            cost_new = (solut.seq[to_1D(0, k, dimen)].C
                        + solut.seq[to_1D(i, j, dimen)].W * cost_concat_1 + solut.seq[to_1D(i, j, dimen)].C
                        + solut.seq[to_1D(k_next, i_prev, dimen)].W * cost_concat_2
                        + solut.seq[to_1D(k_next, i_prev, dimen)].C
                        + solut.seq[to_1D(j_next, dimen, dimen)].W * cost_concat_3
                        + solut.seq[to_1D(j_next, dimen, dimen)].C)

            if cost_new < cost_best:
                cost_best = cost_new - DBL_EPSILON
                I, J, POS = i, j, k

        for k in range(i + opt, dimen):
            k_next = k + 1

            cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]]
            cost_concat_2 = cost_concat_1 + solut.seq[to_1D(j_next, k, dimen)].T + data.cost[solut.s[k]][solut.s[i]]
            cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i, j, dimen)].T + data.cost[solut.s[j]][solut.s[k_next]]

            cost_new = (solut.seq[to_1D(0, i_prev, dimen)].C
                        + solut.seq[to_1D(j_next, k, dimen)].W * cost_concat_1 + solut.seq[to_1D(j_next, k, dimen)].C
                        + solut.seq[to_1D(i, j, dimen)].W * cost_concat_2 + solut.seq[to_1D(i, j, dimen)].C
                        + solut.seq[to_1D(k_next, dimen, dimen)].W * cost_concat_3
                        + solut.seq[to_1D(k_next, dimen, dimen)].C)

            if cost_new < cost_best:
                cost_best = cost_new - DBL_EPSILON
                I, J, POS = i, j, k

        i += 1
        j += 1

    if cost_best < solut.cost - DBL_EPSILON:
        reinsert(solut.s, I, J, POS + 1)
        update_subseq_info_matrix(solut, data, I if I < POS + 1 else POS + 1)
        return True

    return False


def rvnd(solut, data):
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

    while (A_start <= B_start <= A_end) or (B_start <= A_start <= B_end):
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


def gils_rvnd(imax, iils, data, verbose):
    solut_partial = Solution(data)
    solut_crnt = Solution(data)
    solut_best = Solution(data)

    for i in range(imax):
        aux = data.rnd[data.rnd_index]
        data.rnd_index += 1

        alpha = r_table(aux)

        if verbose:
            print(f"[+] Search {i + 1}")
            print("\t[+] Constructing..")

        solut_crnt.s = construct(alpha, data)
        update_subseq_info_matrix(solut_crnt, data)

        solution_copy(solut_crnt, solut_partial)
        if verbose:
            print("\t[+] Looking for the best Neighbor..")
            print(f"\t    Construction Cost: {solut_partial.cost:.3f}")

        iterILS = 0
        while iterILS < iils:
            rvnd(solut_crnt, data)
            if solut_crnt.cost < solut_partial.cost - DBL_EPSILON:
                solution_copy(solut_crnt, solut_partial)
                iterILS = 0

            solut_crnt.s = perturb(solut_partial, data)
            update_subseq_info_matrix(solut_crnt, data)
            iterILS += 1

        if solut_partial.cost < solut_best.cost - DBL_EPSILON:
            solution_copy(solut_partial, solut_best)

        if verbose:
            print(f"\tCurrent best cost: {solut_best.cost}")
            print("SOLUCAO: " + " ".join(str(v) for v in solut_best.s))

    print(f"COST: {solut_best.cost:.2f}")


def main():
    Imax = 10
    verbose = False

    for arg in sys.argv[1:]:
        if arg in ("-v", "--verbose"):
            verbose = True

    dimen, cost, rnd = load_data()
    data = TData(cost, dimen, rnd)

    random.seed(time.process_time())

    Iils = data.dimen if data.dimen < 100 else 100

    t1 = time.perf_counter()
    start = time.process_time()
    gils_rvnd(Imax, Iils, data, verbose)
    cpu_time = time.process_time() - start
    t2 = time.perf_counter()
    duration = (t2 - t1) * 1000

    res = duration / 10e2
    print(f"TIME: {cpu_time}")
    print(f"wall clock (s): {res}")

    print(f"Tamanho RND {len(rnd)}")

    sys.exit(0)


if __name__ == "__main__":
    main()
