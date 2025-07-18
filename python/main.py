#! /usr/bin/pypy3

import time
from read import *
import math
import sys

class tRnd:
    def __init__(self, rnd):
        self.rnd = rnd
        self.rnd_index = 0

class tSubseq:
    def __init__(self, T=0.0, W=0.0, C=0.0):
        self.T = T
        self.C = C
        self.W = W


EPSILON = 1e-13

SWAP        = 0
TWO_OPT     = 1
REINSERTION = 2
OR_OPT_2    = 3
OR_OPT_3    = 4

n, m, rnd = get_instance_info()

W = 0
T = 1
C = 2

IT = 0

t_reinsertion = 0
t_or_opt2 = 0
t_or_opt3 = 0
t_swap = 0
t_two_opt = 0
t_construction = 0
t_perturb = 0
t_seq = 0

improv_flag = None

def subseq_info_fill(n):
    matrix = []

    for i in range(n+1):
        matrix.append([])
        for j in range(n+1):
            matrix[i].append(tSubseq(0.0, 0.0, 0.0))
    return matrix

def quicksort(arr, left, right, r):
    if left < right:
        pivot = partition(arr, left, right, r)
        quicksort(arr, left, pivot - 1, r)
        quicksort(arr, pivot + 1, right, r)

def partition(arr, left, right, r):
    pivot = arr[right]
    i = left - 1
    for j in range(left, right):
        if m[r][arr[j]] < m[r][pivot]:
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1

def sort(arr, r):
    quicksort(arr, 0, len(arr) - 1, r)

def construction(alpha, rnd):
    s = [0]
    c_list = list(range(1, n))

    r = 0
    while len(c_list) > 0:

        sort(c_list, r)

        index = rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1
        c = c_list[index]

        s.append(c)
        r = c
        c_list.remove(c)

    s.append(0)

    return s

def update_subseq_info_matrix(sol):
    i = 0
    d = n + 1
    while i < d:
        k = 1 - i - int(not i)

        subseq[i][i].T = 0.0
        subseq[i][i].C = 0.0
        subseq[i][i].W = int(not (i == 0))

        j = i + 1
        while j < d:
            a = j - 1

            subseq[i][j].T = m[sol[a]][sol[j]] + subseq[i][a].T
            subseq[i][j].C = subseq[i][j].T  + subseq[i][a].C
            subseq[i][j].W = j + k

            j += 1

        i += 1

def swap(s, i, j):
    s[i], s[j] = s[j], s[i]

def reverse(s, i, j):
    s[i:j+1] = s[i:j+1][::-1]

def reinsert(s, i, j, pos):
    if i < pos:
        s[pos:pos] = s[i:j+1]
        s[:] = s[:i] + s[j+1:]
    else:
        sub = s[i:j+1]
        s[:] = s[:i] + s[j+1:]
        s[pos:pos] = sub


def search_swap(s):
    global subseq
    cost_best = float('inf')
    I = None
    J = None

    for i in range(1, (n-1) ):
        i_prev = i - 1
        i_next = i + 1

        cost_concat_1 = subseq[0][i_prev].T  + m[s[i_prev]][s[i_next]]
        cost_concat_2 = cost_concat_1 + subseq[i][i_next].T + m[s[i]][s[i_next+1]]

        cost = subseq[0][i_prev].C                                            + \
                subseq[i][i_next].W   * (cost_concat_1) + m[s[i_next]][s[i]]  + \
                subseq[i_next+1][n].W * (cost_concat_2) + subseq[i_next+1][n].C

        if cost < cost_best:
            cost_best = cost - EPSILON
            I = i
            J = i_next

        for j in range(i_next+1, (n)):
            j_next = j + 1
            j_prev = j - 1

            cost_concat_1 = subseq[0][i_prev].T + m[s[i_prev]][s[j]]
            cost_concat_2 = cost_concat_1 + m[s[j]][s[i_next]]
            cost_concat_3 = cost_concat_2 + subseq[i_next][j_prev].T + m[s[j_prev]][s[i]]
            cost_concat_4 = cost_concat_3 + m[s[i]][s[j_next]]

            cost = subseq[0][i_prev].C                                                + \
                    cost_concat_1                                                   + \
                    subseq[i_next][j_prev].W * cost_concat_2 + subseq[i_next][j_prev].C + \
                    cost_concat_3                                                   + \
                    subseq[j_next][n].W * cost_concat_4 + subseq[j_next][n].C

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j

    if cost_best < subseq[0][n].C  - EPSILON:
        swap(s, I, J)
        update_subseq_info_matrix(s)
        global improv_flag
        improv_flag = True

def search_two_opt(s):
    global subseq
    cost_best = float('inf')


    for i in range(1, n-1):
        i_prev = i - 1
        
        reverse_cost = subseq[i][i+1].T 
        for j in range(i+2, n):
            j_next = j + 1

            reverse_cost += m[s[j-1]][s[j]] * (subseq[i][j].W -1)

            cost_concat_1 = subseq[0][i_prev].T  + m[s[j]][s[i_prev]]
            cost_concat_2 = cost_concat_1 + subseq[i][j].T + m[s[j_next]][s[i]]

            cost = subseq[0][i_prev].C                                        + \
                    subseq[i][j].W * cost_concat_1 + reverse_cost             + \
                    subseq[j_next][n].W  * (cost_concat_2) + subseq[j_next][n].C 

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j

    if cost_best < subseq[0][n].C  - EPSILON:
        reverse(s, I, J)
        update_subseq_info_matrix(s)
        global improv_flag
        improv_flag = True


def search_reinsertion(s, OPT):
    global subseq
    cost_best = float('inf')
    I = None
    J = None
    POS = None
    opt = OPT - 1 
    MAX = n - opt 

    for i in range(1, MAX + 1):
        j = opt + i - 1

        j_next = j+1
        i_prev = i-1

        for k in range(0, i_prev):
            k_next = k+1

            cost_concat_1 = subseq[0][k].T  + m[s[k]][s[i]]
            cost_concat_2 = cost_concat_1 + subseq[i][j].T  + m[s[j]][s[k_next]]
            cost_concat_3 = cost_concat_2 + subseq[k_next][i_prev].T  + m[s[i_prev]][s[j_next]]

            cost = subseq[0][k].C                                                          + \
                    subseq[i][j].W             * cost_concat_1 + subseq[i][j].C               + \
                    subseq[k_next][i_prev].W   * cost_concat_2 + subseq[k_next][i_prev].C     + \
                    subseq[j_next][n].W        * cost_concat_3 + subseq[j_next][n].C 

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j
                POS = k

        for k in range(i+opt, n):
            k_next = k+1

            cost_concat_1 = subseq[0][i_prev].T  + m[s[i_prev]][s[j_next]]
            cost_concat_2 = cost_concat_1 + subseq[j_next][k].T  + m[s[k]][s[i]]
            cost_concat_3 = cost_concat_2 + subseq[i][j].T  + m[s[j]][s[k_next]]

            cost = subseq[0][i_prev].C                                         + \
                    subseq[j_next][k].W    * cost_concat_1 + subseq[j_next][k].C  + \
                    subseq[i][j].W         * cost_concat_2 + subseq[i][j].C       + \
                    subseq[k_next][n].W    * cost_concat_3 + subseq[k_next][n].C 

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j
                POS = k

    if cost_best < subseq[0][n].C  - EPSILON:
        reinsert(s, I, J, POS+1)
        update_subseq_info_matrix(s)
        global improv_flag
        improv_flag = True

def RVND(s, rnd):

    global subseq
    global improv_flag
    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]

    while len(neighbd_list) > 0:
        i = rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1

        neighbd = neighbd_list[i]

        improv_flag = False

        if neighbd == SWAP:
            search_swap(s)
        elif neighbd == TWO_OPT:
            search_two_opt(s)
        elif neighbd == REINSERTION:
            search_reinsertion(s, REINSERTION)
        elif neighbd == OR_OPT_2:
            search_reinsertion(s, OR_OPT_2)
        elif neighbd == OR_OPT_3:
            search_reinsertion(s, OR_OPT_3)

        if improv_flag == True:
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]
        else:
            neighbd_list.pop(i)

    return s

def perturb(sl, rnd):
    s = sl[:]
    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end):

        A_start = rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1
        A_end = A_start + rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1

        B_start = rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1
        B_end = B_start + rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1

    if A_start < B_start:
        reinsert(s, B_start, B_end-1, A_end)
        reinsert(s, A_start, A_end-1, B_end )
    else:
        reinsert(s, A_start, A_end-1, B_end)
        reinsert(s, B_start, B_end-1, A_end )

    return s

def GILS_RVND(Imax, Iils, R, rnd):

    global subseq

    cost_best = float('inf')
    s_best = []


    #  TODO:  tornar subseq global
    subseq = subseq_info_fill(n)

    for i in range(Imax):
        index = rnd.rnd[rnd.rnd_index]
        rnd.rnd_index += 1

        alpha = R[index]

        print("[+] Local Search {}".format(i+1))
        print("\t[+] Constructing Inital Solution..")
        s = construction(alpha, rnd)
        update_subseq_info_matrix(s)
        sl = s[:]
        rvnd_cost_best = subseq[0][n].C  - EPSILON

        print(rvnd_cost_best)

        print("\t[+] Looking for the best Neighbor..")
        iterILS = 0
        while iterILS < Iils:
            RVND(s, rnd)
            rvnd_cost_crnt  = subseq[0][n].C  - EPSILON
            if rvnd_cost_crnt < rvnd_cost_best:
                rvnd_cost_best = rvnd_cost_crnt
                sl = s[:]
                iterILS = 0

            s = perturb(sl, rnd)
            update_subseq_info_matrix(s)
            iterILS += 1

        update_subseq_info_matrix(sl)
        sl_cost = subseq[0][n].C  - EPSILON

        if sl_cost < cost_best:
            s_best = sl[:]
            cost_best = sl_cost

        print("\tCurrent best solution cost: {}".format(cost_best))

    print("COST: {}".format (cost_best))
    print("SOLUTION: {}".format( s_best))


def main(rnd):

    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]

    Imax = 10
    Iils = min(n, 100)

    GILS_RVND(Imax, Iils, R, tRnd(rnd))

start = time.time()
main(rnd)
print("TIME: %s " % (time.time() - start))
