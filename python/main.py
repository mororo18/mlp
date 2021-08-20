#! /usr/bin/pypy

import time
from read import *
from random import randint
import math

EPSILON = 1e-13

SWAP        = 0
TWO_OPT     = 1
REINSERTION = 2
OR_OPT_2    = 3
OR_OPT_3    = 4

n, m = get_instance_info()

W = "W"
T = "T"
C = "C"

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
    matrix = {
            W:[],
            T:[],
            C:[]
            }

    for i in range(n+1):
        matrix[W].append([])
        matrix[T].append([])
        matrix[C].append([])
        for j in range(n+1):
            matrix[W][i].append(0.0)
            matrix[T][i].append(0.0)
            matrix[C][i].append(0.0)

    return matrix

def construction(alpha):
    s = [0]
    c_list = list(range(1, n))

    r = 0
    while len(c_list) > 0:

        i = int(len(c_list)*alpha) + 1

        c_list = sorted(c_list, key = lambda n : m[n][r], reverse=False)

        c = c_list[randint(0, i-1)]
        s.append(c)
        r = c
        c_list.remove(c)

    s.append(0)

    return s

def subseq_info_load(sol, seq):
    global t_seq
    t_seq -= time.time()
    i = 0
    d = n + 1
    while i < d:
        k = 1 - i - int(not i)

        seq[T][i][i] = 0.0
        seq[C][i][i] = 0.0
        seq[W][i][i] = int(not (i == 0))

        j = i + 1
        while j < d:
            a = j - 1

            seq[T][i][j] = m[sol[a]][sol[j]] + seq[T][i][a]
            seq[C][i][j] = seq[T][i][j] + seq[C][i][a]
            seq[W][i][j] = j + k

            """
            seq[T][j][i] = seq[T][i][j]
            seq[C][j][i] = seq[C][i][j]
            seq[W][j][i] = seq[W][i][j]
            """

            j += 1

        i += 1

    t_seq += time.time()

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


def search_swap(s, seq):
    cost_best = float('inf')
    I = None
    J = None

    for i in range(1, (n-1) ):
        i_prev = i - 1
        i_next = i + 1

        cost_concat_1 = seq[T][0][i_prev] + m[s[i_prev]][s[i_next]]
        cost_concat_2 = cost_concat_1 + seq[T][i][i_next] + m[s[i]][s[i_next+1]]

        cost = seq[C][0][i_prev]                                            + \
                seq[W][i][i_next]   * (cost_concat_1) + m[s[i_next]][s[i]]  + \
                seq[W][i_next+1][n] * (cost_concat_2) + seq[C][i_next+1][n]

        if cost < cost_best:
            cost_best = cost - EPSILON
            I = i
            J = i_next

        for j in range(i_next+1, (n)):
            j_next = j + 1
            j_prev = j - 1

            cost_concat_1 = seq[T][0][i_prev] + m[s[i_prev]][s[j]]
            cost_concat_2 = cost_concat_1 + m[s[j]][s[i_next]]
            cost_concat_3 = cost_concat_2 + seq[T][i_next][j_prev] + m[s[j_prev]][s[i]]
            cost_concat_4 = cost_concat_3 + m[s[i]][s[j_next]]

            """
            cost = 1st subseq                       +
                    concat 2nd subseq (sigle node)  + 
                    concat 3rd subseq               +
                    concat 4th subseq (single node) +
                    concat 5th subseq
            """

            cost = seq[C][0][i_prev]                                                + \
                    cost_concat_1                                                   + \
                    seq[W][i_next][j_prev] * cost_concat_2 + seq[C][i_next][j_prev] + \
                    cost_concat_3                                                   + \
                    seq[W][j_next][n] * cost_concat_4 + seq[C][j_next][n]

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j

    if cost_best < seq[C][0][n] - EPSILON:
        #print(cost_best, I, J)
        swap(s, I, J)
        subseq_info_load(s, seq)
        #print(seq[C][0][n])
        global improv_flag
        improv_flag = True

def search_two_opt(s, seq):
    cost_best = float('inf')

    for i in range(1, n-1):
        i_prev = i - 1
        
        reverse_cost = seq[T][i][i+1]
        for j in range(i+2, n):
            j_next = j + 1

            reverse_cost += m[s[j-1]][s[j]] * (seq[W][i][j]-1)

            cost_concat_1 = seq[T][0][i_prev] + m[s[j]][s[i_prev]]
            cost_concat_2 = cost_concat_1 + seq[T][i][j] + m[s[j_next]][s[i]]

            """
            cost = 1st subseq                           +
                    concat 2nd subseq (reversed seq)    + 
                    concat 3rd subseq               
            """

            cost = seq[C][0][i_prev]                                        + \
                    seq[W][i][j] * cost_concat_1 + reverse_cost             + \
                    seq[W][j_next][n] * (cost_concat_2) + seq[C][j_next][n]

            #print(cost, i, "        ", j)
            #print(seq[C][0][i_prev], seq[W][i][j] * cost_concat_1,seq[W][j_next][n] * (cost_concat_2), seq[C][j_next][n] )
            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j

    if cost_best < seq[C][0][n] - EPSILON:
        #print(cost_best)
        #print(I, J)
        reverse(s, I, J)
        subseq_info_load(s, seq)
        global improv_flag
        improv_flag = True


def search_reinsertion(s, seq, OPT):
    cost_best = float('inf')
    I = None
    J = None
    POS = None
    opt = OPT - 1 
    MAX = n - opt 

    '''
    #branchless approach
    BEST = 1
    J_arr = [None, None]
    I_arr = [None, None]
    POS_arr = [None, None]
    COST_arr=[None, cost_best]
    '''

    for i in range(1, MAX + 1):
        j = opt + i - 1

        j_next = j+1
        i_prev = i-1

        for k in range(0, i_prev):
            k_next = k+1

            cost_concat_1 = seq[T][0][k] + m[s[k]][s[i]]
            cost_concat_2 = cost_concat_1 + seq[T][i][j] + m[s[j]][s[k_next]]
            cost_concat_3 = cost_concat_2 + seq[T][k_next][i_prev] + m[s[i_prev]][s[j_next]]

            """
            cost = 1st subseq                           +
                    concat 2nd subseq (reinserted seq)  + 
                    concat 3rd subseq                   +
                    concat 4th subseq  
            """

            cost = seq[C][0][k]                                                         + \
                    seq[W][i][j]            * cost_concat_1 + seq[C][i][j]              + \
                    seq[W][k_next][i_prev]  * cost_concat_2 + seq[C][k_next][i_prev]    + \
                    seq[W][j_next][n]       * cost_concat_3 + seq[C][j_next][n]

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j
                POS = k
            '''
            #branchless approach
            if_state = int(cost < COST_arr[BEST])
            COST_arr[if_state] = cost - EPSILON
            I_arr[if_state] = i
            J_arr[if_state] = j
            POS_arr[if_state] = k
            '''


        for k in range(i+opt, MAX - 1):
            k_next = k+1

            cost_concat_1 = seq[T][0][i_prev] + m[s[i_prev]][s[j_next]]
            cost_concat_2 = cost_concat_1 + seq[T][j_next][k] + m[s[k]][s[i]]
            cost_concat_3 = cost_concat_2 + seq[T][i][j] + m[s[j]][s[k_next]]

            """
            cost = 1st subseq                           +
                    concat 2nd subseq                   + 
                    concat 3rd subseq  (reinserted seq) +
                    concat 4th subseq  
            """

            cost = seq[C][0][i_prev]                                        + \
                    seq[W][j_next][k]   * cost_concat_1 + seq[C][j_next][k] + \
                    seq[W][i][j]        * cost_concat_2 + seq[C][i][j]      + \
                    seq[W][k_next][n]   * cost_concat_3 + seq[C][k_next][n]

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j
                POS = k
            '''
            #branchless approach
            if_state = int(cost < COST_arr[BEST])
            COST_arr[if_state] = cost - EPSILON
            I_arr[if_state] = i
            J_arr[if_state] = j
            POS_arr[if_state] = k
            '''

    #if COST_arr[BEST] < seq[C][0][n] - EPSILON:
    if cost_best < seq[C][0][n] - EPSILON:
        #reinsert(s, I_arr[BEST], J_arr[BEST], POS_arr[BEST]+1)
        reinsert(s, I, J, POS+1)
        subseq_info_load(s, seq)
        global improv_flag
        improv_flag = True

def RVND(s, subseq):


    global t_reinsertion
    global t_or_opt2 
    global t_or_opt3 
    global t_swap 
    global t_two_opt 

    global improv_flag
    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]

    '''
        s = []
        
        s.append(0)
        i = n-1
        while i > 0:
            s.append(i)
            i -= 1
        s.append(0)
        print(s)
        subseq_info_load(s, subseq)

        print(subseq[C][0][n])
    '''

    while len(neighbd_list) > 0:
        global IT
        IT += 1
        i = randint(0, len(neighbd_list)-1)
        neighbd = neighbd_list[i]

        improv_flag = False

        if neighbd == SWAP:
            t_swap -= time.time()
            search_swap(s, subseq)
            t_swap += time.time()
        elif neighbd == TWO_OPT:
            t_two_opt -= time.time()
            search_two_opt(s, subseq)
            t_two_opt += time.time()
        elif neighbd == REINSERTION:
            t_reinsertion -= time.time()
            search_reinsertion(s, subseq, REINSERTION)
            t_reinsertion += time.time()
        elif neighbd == OR_OPT_2:
            t_or_opt2 -= time.time()
            search_reinsertion(s, subseq, OR_OPT_2)
            t_or_opt2 += time.time()
        elif neighbd == OR_OPT_3:
            t_or_opt3 -= time.time()
            search_reinsertion(s, subseq, OR_OPT_3)
            t_or_opt3 += time.time()


        if improv_flag == True:
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]
        else:
            neighbd_list.pop(i)

    return s

def perturb(sl):
    s = sl[:]
    #s = sl.copy()
    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

    size_max = math.floor(len(s)/10) if math.floor(len(s)/10) >= 2 else 2
    size_min = 2

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end):
        A_start = randint(1, len(s) - 1 - size_max)
        A_end = A_start + randint(size_min, size_max)

        B_start = randint(1, len(s) - 1 - size_max)
        B_end = B_start + randint(size_min, size_max)

    if A_start < B_start:
        reinsert(s, B_start, B_end-1, A_end)
        reinsert(s, A_start, A_end-1, B_end )
    else:
        reinsert(s, A_start, A_end-1, B_end)
        reinsert(s, B_start, B_end-1, A_end )


    return s

def GILS_RVND(Imax, Iils, R):

    cost_best = float('inf')
    s_best = []

    global t_construction
    global t_perturb

    subseq = subseq_info_fill(n)

    for i in range(Imax):
        alpha = R[randint(0, len(R)-1)]

        print("[+] Local Search {}".format(i+1))
        print("\t[+] Constructing Inital Solution..")
        t_construction -= time.time()
        s = construction(alpha)
        t_construction += time.time()
        subseq_info_load(s, subseq)
        sl = s[:]
        #sl = s.copy()
        rvnd_cost_best = subseq[C][0][n] - EPSILON

        print("\t[+] Looking for the best Neighbor..")
        iterILS = 0
        while iterILS < Iils:
            RVND(s, subseq)
            rvnd_cost_crnt  = subseq[C][0][n] - EPSILON
            if rvnd_cost_crnt < rvnd_cost_best:
                rvnd_cost_best = rvnd_cost_crnt
                sl = s[:]
                #sl = s.copy()
                iterILS = 0

            t_perturb -= time.time()
            s = perturb(sl)
            t_perturb += time.time()
            subseq_info_load(s, subseq)
            iterILS += 1

        subseq_info_load(sl, subseq)
        sl_cost = subseq[C][0][n] - EPSILON

        if sl_cost < cost_best:
            s_best = sl
            cost_best = sl_cost

        print("\tCurrent best solution cost: {}".format(cost_best))

    print("COST: {}".format (cost_best))
    print("SOLUTION: {}".format( s_best))
    print("Total Iterations RVND {}".format(IT))




def main():

    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]

    Imax = 10
    Iils = min(n, 100)

    GILS_RVND(Imax, Iils, R)

start = time.time()
main()
print("TIME %s " % (time.time() - start))
print("SWAP %s" % t_swap)
print("Reinsert %s" % t_reinsertion)
print("or_opt2 %s" % t_or_opt2)
print("or_opt3 %s" % t_or_opt3)
print("two_opt %s" % t_two_opt)
print("subseq_load %s" % t_seq)
