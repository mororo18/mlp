#! /usr/bin/pypy3

import time
from read import *
import math
import sys
#import matplotlib.pyplot as plt 
from typing import List, NoReturn

Matrix = List[List[float]]

class tInfo:
    IT : int = 0
    dimen : int
    cost : Matrix

    T : int = 0
    C : int = 1
    W : int = 2

    SWAP : int = 0
    REINSERTION : int = 1
    OR_OPT_2 : int = 2
    OR_OPT_3 : int = 3
    TWO_OPT : int = 4

    EPSILON : float = 1e-13

    rnd : List[int]
    rnd_index : int = 0

    def __init__(self, dimension : int, cost : Matrix, rnd : List[int]):
        self.dimen = dimension
        self.cost = cost
        self.rnd = rnd

class tSolution:
    s : List[int]
    seq : List[List[List[float]]]
    cost : float

    def __init__(self, s : List[int], seq : List[List[List[float]]], cost : float):
        self.s = s
        self.seq = seq
        self.cost = cost


def subseq_fill(info : tInfo) -> List[List[List[float]]]:
    arr = []

    for i in range(((info.dimen+1)**2)*3):
        arr.append(0.0)
       #arr.append([])
       #for j in range(info.dimen+1):
       #    arr[i].append([0.0, 0.0, 0.0])

    return arr

def quicksort(arr, left, right, info, r):
    if left < right:
        pivot = partition(arr, left, right, info, r)
        quicksort(arr, left, pivot - 1, info, r)
        quicksort(arr, pivot + 1, right, info, r)

def partition(arr, left, right, info, r):
    pivot = arr[right]
    i = left - 1
    for j in range(left, right):
        if info.cost[r][arr[j]] < info.cost[r][pivot]:
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1

def sort(arr, r, info):
    quicksort(arr, 0, len(arr) - 1, info, r)

def construction(alpha : float, info : tInfo) -> List[int]:
    s : List[int] = [0]
    c_list : List[int] = list(range(1, info.dimen))

    r = 0
    while len(c_list) > 0:
        sort(c_list, r, info)

        r_value = info.rnd[info.rnd_index]
        info.rnd_index += 1

        c = c_list[r_value]

        s.append(c)
        r = c
        c_list.remove(c)

    s.append(0)

    return s

def to_1D(i, j, k, dimen):
    return 3*(i*dimen + j) + k

def subseq_load(solut : tSolution, info : tInfo) -> NoReturn:
    for i in range(0, info.dimen+1):
        k : int = 1 - i - int(not i)

        solut.seq[to_1D(i, i, info.T, info.dimen+1)] = 0.0
        solut.seq[to_1D(i, i, info.C, info.dimen+1)] = 0.0
        solut.seq[to_1D(i, i, info.W, info.dimen+1)] = float(int(not (i == 0)))

        for j in range(i+1, info.dimen+1):
        #while j < d:
            #print(i, j)
            j_prev : int = j - 1

            solut.seq[to_1D(i, j, info.T, info.dimen+1)] = info.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[to_1D(i, j_prev, info.T, info.dimen+1)]
            solut.seq[to_1D(i, j, info.C, info.dimen+1)] = solut.seq[to_1D(i, j, info.T, info.dimen+1)] + solut.seq[to_1D(i, j_prev, info.C, info.dimen+1)]
            solut.seq[to_1D(i, j, info.W, info.dimen+1)] = float(j + k)

        #exit(0)

    solut.cost = solut.seq[to_1D(0, info.dimen, info.C, info.dimen+1)] - info.EPSILON

def swap(s : List[int], i : int, j : int) -> NoReturn:
    s[i], s[j] = s[j], s[i]

def reverse(s : List[int], i : int, j : int) -> NoReturn:
    s[i:j+1] = s[i:j+1][::-1]

def reinsert(s : List[int], i : int, j : int, pos : int) -> NoReturn:
    if i < pos:
        s[pos:pos] = s[i:j+1]
        s[:] = s[:i] + s[j+1:]
    else:
        sub = s[i:j+1]
        s[:] = s[:i] + s[j+1:]
        s[pos:pos] = sub


def search_swap(solut : tSolution, info : tInfo) -> bool:
    cost_new : float
    cost_concat_1 : float
    cost_concat_2 : float
    cost_concat_3 : float
    cost_concat_4 : float

    cost_best : float = float('inf')
    I : int = -1
    J : int = -1

    for i in range(1, (info.dimen-1)):
        i_prev : int = i - 1
        i_next : int = i + 1

        cost_concat_1 = solut.seq[to_1D(0, i_prev, info.T, info.dimen+1)] + info.cost[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, i_next, info.T, info.dimen+1)] + info.cost[solut.s[i]][solut.s[i_next+1]]

        cost_new = solut.seq[to_1D(0, i_prev, info.C, info.dimen+1)]                                            + \
                solut.seq[to_1D(i, i_next, info.W, info.dimen+1)]   * (cost_concat_1) + info.cost[solut.s[i_next]][solut.s[i]]  + \
                solut.seq[to_1D(i_next+1, info.dimen, info.W, info.dimen+1)] * (cost_concat_2) + solut.seq[to_1D(i_next+1, info.dimen, info.C, info.dimen+1)]

        if cost_new < cost_best:
            cost_best = cost_new - info.EPSILON
            I = i
            J = i_next
            #print(cost_best, I, J)

        for j in range(i_next+1, (info.dimen)):
            j_next : int = j + 1
            j_prev : int = j - 1

            cost_concat_1 = solut.seq[to_1D(0, i_prev, info.T, info.dimen+1)] + info.cost[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + info.cost[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i_next, j_prev, info.T, info.dimen+1)] + info.cost[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3 + info.cost[solut.s[i]][solut.s[j_next]]

            cost_new = solut.seq[to_1D(0, i_prev, info.C, info.dimen+1)]                                                + \
                    cost_concat_1                                                   + \
                    solut.seq[to_1D(i_next, j_prev, info.W, info.dimen+1)] * cost_concat_2 + solut.seq[to_1D(i_next, j_prev, info.C, info.dimen+1)] + \
                    cost_concat_3                                                   + \
                    solut.seq[to_1D(j_next, info.dimen, info.W, info.dimen+1)] * cost_concat_4 + solut.seq[to_1D(j_next, info.dimen, info.C, info.dimen+1)]

            if cost_new < cost_best:
                cost_best = cost_new - info.EPSILON
                I = i
                J = j
                #print(cost_best, I, J)

    #print(cost_best, solut.cost, I, J)
    if cost_best < solut.cost - info.EPSILON:
        swap(solut.s, I, J)
        subseq_load(solut, info)
        return True
    #print(seq[C][0][n])
        #if cost_best != seq[C][0][n]:
        #   print("ERRR")
       #global improv_flag
       #global improv_swap
       #print("\nswap")
       #print(cost_best)
       #print(solut.cost)

    return False

def search_two_opt(solut : tSolution, info : tInfo) -> bool:
    cost_new : float
    cost_concat_1 : float
    cost_concat_2 : float
    cost_best = float('inf')
    I : int
    J : int

    for i in range(1, info.dimen-1):
        i_prev : int = i - 1
        
        reverse_cost : float = solut.seq[to_1D(i, i+1, info.T, info.dimen+1)]
        for j in range(i+2, info.dimen):
            j_next : int = j + 1

            reverse_cost += info.cost[solut.s[j-1]][solut.s[j]] * (solut.seq[to_1D(i, j, info.W, info.dimen+1)]-1)

            cost_concat_1 = solut.seq[to_1D(0, i_prev, info.T, info.dimen+1)] + info.cost[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, info.T, info.dimen+1)] + info.cost[solut.s[j_next]][solut.s[i]]

            cost_new = solut.seq[to_1D(0, i_prev, info.C, info.dimen+1)]                                        + \
                    solut.seq[to_1D(i, j, info.W, info.dimen+1)] * cost_concat_1 + reverse_cost             + \
                    solut.seq[to_1D(j_next, info.dimen, info.W, info.dimen+1)] * (cost_concat_2) + solut.seq[to_1D(j_next, info.dimen, info.C, info.dimen+1)]

            #print(cost, i, "        ", j)
            #print(seq[to_1D(C, 0, i_prev, info.dimen+1)], seq[W][i][j] * cost_concat_1,seq[W][j_next][n] * (cost_concat_2), seq[C][j_next][n] )
            if cost_new < cost_best:
                cost_best = cost_new - info.EPSILON
                I = i
                J = j

            #if i == 2 and j == 10:
                #print ("opaa")

    if cost_best < solut.cost - info.EPSILON:
        #print(cost_best)
        #print(I, J)
        reverse(solut.s, I, J)
        subseq_load(solut, info)
#       global improv_flag
#       global improv_two_opt
#       improv_flag = True
#       improv_two_opt += 1
       #print("\ntwo_opt")
       #print(cost_best)
       #print(solut.cost)
        return True

    return False

"""
Movimento avalido me maneira incorreta (nesse caso). Os que presenciei aparentavam realizar a avaliacao de maneira correta
[0, 7, 8, 10, 12, 6, 5, 4, 1, 2, 3, 11, 13, 9, 0]
[0, 7, 8, 10, 12, 6, 11, 5, 4, 1, 2, 3, 13, 9, 0]

reinsertion 1
15752.0
24694.0
11 11 5
"""

# TODO (EUEU) testar acessando as infos das estruturas pelo 1o indice do array3D. Ex: seq[C][...

#def search_reinsertion(solut , info , opt , seq ) :
def search_reinsertion(solut : tSolution, info : tInfo, opt : int, seq : List[List[List[float]]]) -> bool:
    """
    cost_concat_1 : float
    cost_concat_2 : float
    cost_concat_3 : float
    cost_new : float
    """

    cost_best : float = float('inf')
    #cost_best : float = float('inf')

    """
    I : int 
    J : int 
    POS : int
    MAX : int = info.dimen - opt 
    """
    I = -1
    J  = -1 
    POS  = -1
    MAX = info.dimen - opt 


    #seq = solut.seq

    for i in range(1, MAX + 1):
        j  = opt + i - 1
        #j : int = opt + i - 1

        j_next  = j+1
        i_prev  = i-1
       #j_next : int = j+1
       #i_prev : int = i-1

        for k in range(0, i_prev):
            k_next  = k+1
            #k_next : int = k+1


            cost_concat_1 =                 seq[to_1D(0, k, info.T, info.dimen+1)] + info.cost[solut.s[k]][solut.s[i]]
            cost_concat_2 = cost_concat_1 + seq[to_1D(i, j, info.T, info.dimen+1)] + info.cost[solut.s[j]][solut.s[k_next]]
            cost_concat_3 = cost_concat_2 + seq[to_1D(k_next, i_prev, info.T, info.dimen+1)] + info.cost[solut.s[i_prev]][solut.s[j_next]]


            cost_new = seq[to_1D(0, k, info.C, info.dimen+1)]                                                                     + \
                    seq[to_1D(i, j, info.W, info.dimen+1)]  * cost_concat_1 + seq[to_1D(i, j, info.C, info.dimen+1)] + \
                    seq[to_1D(k_next, i_prev, info.W, info.dimen+1)]  * cost_concat_2 + seq[to_1D(k_next, i_prev, info.C, info.dimen+1)] + \
                    seq[to_1D(j_next, info.dimen, info.W, info.dimen+1)]  * cost_concat_3 + seq[to_1D(j_next, info.dimen, info.C, info.dimen+1)]

            if cost_new < cost_best:
                cost_best = cost_new - info.EPSILON
                I = i
                J = j
                POS = k
            #if i == 11 and k == 5:
                #print ("opa")


        for k in range(i+opt, info.dimen):
            k_next = k+1
            #k_next : int = k+1

            cost_concat_1 =                 seq[to_1D(0, i_prev, info.T, info.dimen+1)] + info.cost[solut.s[i_prev]][solut.s[j_next]]
            cost_concat_2 = cost_concat_1 + seq[to_1D(j_next, k, info.T, info.dimen+1)] + info.cost[solut.s[k]]     [solut.s[i]]
            cost_concat_3 = cost_concat_2 + seq[to_1D(i, j, info.T, info.dimen+1)] + info.cost[solut.s[j]]     [solut.s[k_next]]


            cost_new = seq[to_1D(0, i_prev, info.C, info.dimen+1)]                                        + \
                    seq[to_1D(j_next, k, info.W, info.dimen+1)]   * cost_concat_1 + seq[to_1D(j_next, k, info.C, info.dimen+1)] + \
                    seq[to_1D(i, j, info.W, info.dimen+1)]        * cost_concat_2 + seq[to_1D(i, j, info.C, info.dimen+1)]      + \
                    seq[to_1D(k_next, info.dimen, info.W, info.dimen+1)]   * cost_concat_3 + seq[to_1D(k_next, info.dimen, info.C, info.dimen+1)]

            if cost_new < cost_best:
                cost_best = cost_new - info.EPSILON
                I = i
                J = j
                POS = k

    #if COST_arr[BEST] < seq[C][0][n] - EPSILON:
    if cost_best < solut.cost - info.EPSILON:
        #reinsert(s, I_arr[BEST], J_arr[BEST], POS_arr[BEST]+1)
        #print(solut.s)
        reinsert(solut.s, I, J, POS+1)
        #print(solut.s)
        subseq_load(solut, info)
        """
        print("\nreinsertion", opt)
        print(cost_best)
        print(solut.cost)
        """

       #if(cost_best != solut.cost):
       #    print(I, J, POS)
       #    exit(0)
        return True

    return False

def RVND(solut : tSolution, info : tInfo) -> NoReturn:

   #global t_reinsertion
   #global t_or_opt2 
   #global t_or_opt3 
   #global t_swap 
   #global t_two_opt 

   #global improv_flag
    #neighbd_list = [ TWO_OPT]
    neighbd_list : List[int] = [info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT_2, info.OR_OPT_3]

    #print(solut.s)

    while len(neighbd_list) > 0:
        i = info.rnd[info.rnd_index]
        #print(info.rnd_index, i )
        #print(neighbd_list)
        info.rnd_index += 1


        neighbd = neighbd_list[i]

        improve : bool = False

        if neighbd == info.SWAP:
            #t_swap -= time.time()
            improve = search_swap(solut, info)
            #t_swap += time.time()
        elif neighbd == info.TWO_OPT:
           #t_two_opt -= time.time()
            improve = search_two_opt(solut, info)
           #t_two_opt += time.time()
        elif neighbd == info.REINSERTION:
           #t_reinsertion -= time.time()
            improve = search_reinsertion(solut, info, info.REINSERTION, solut.seq)
           #t_reinsertion += time.time()
        elif neighbd == info.OR_OPT_2:
           #t_or_opt2 -= time.time()
            improve = search_reinsertion(solut, info, info.OR_OPT_2, solut.seq)
           #t_or_opt2 += time.time()
        elif neighbd == info.OR_OPT_3:
           #t_or_opt3 -= time.time()
            improve = search_reinsertion(solut, info, info.OR_OPT_3, solut.seq)
           #t_or_opt3 += time.time()

        #print(improve)

        if improve == True:
            neighbd_list = [info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT_2, info.OR_OPT_3]
            
        else:
            neighbd_list.pop(i)

def perturb(sl : List[int], info : tInfo) -> List[int]:
    s = sl.copy()

    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

    size_max = math.floor(len(s)/10) if math.floor(len(s)/10) >= 2 else 2
    size_min = 2

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end):

        A_start = info.rnd[info.rnd_index]
        info.rnd_index += 1
        A_end = A_start + info.rnd[info.rnd_index]
        info.rnd_index += 1

        B_start = info.rnd[info.rnd_index]
        info.rnd_index += 1
        B_end = B_start + info.rnd[info.rnd_index]
        info.rnd_index += 1

    if A_start < B_start:
        reinsert(s, B_start, B_end-1, A_end)
        reinsert(s, A_start, A_end-1, B_end )
    else:
        reinsert(s, A_start, A_end-1, B_end)
        reinsert(s, B_start, B_end-1, A_end )

    return s

def GILS_RVND(Imax : int, Iils : int, R : List[float], info : tInfo) -> NoReturn:

    solut_partial : tSolution = tSolution([0 for i in range(info.dimen+1)], subseq_fill(info), 0.0)
    solut_crnt : tSolution = tSolution([0 for i in range(info.dimen+1)], subseq_fill(info), 0.0)
    solut_best : tSolution = tSolution([0 for i in range(info.dimen+1)], subseq_fill(info), float('inf'))

    for i in range(Imax):
        r_value = info.rnd[info.rnd_index]
        print(r_value)
        info.rnd_index += 1

        alpha : float = R[r_value]

        print("[+] Local Search {}".format(i+1))
        solut_crnt.s = construction(alpha, info)
        subseq_load(solut_crnt, info)

        solut_partial.s = solut_crnt.s.copy()
        solut_partial.cost = solut_crnt.cost
        print("\t[+] Constructing Inital Solution..", solut_crnt.cost)

        print("\t[+] Looking for the best Neighbor..")
        iterILS = 0
        while iterILS < Iils:
            #print("ILS")
            RVND(solut_crnt, info)

            if solut_crnt.cost < solut_partial.cost:
                solut_partial.s = (solut_crnt.s).copy()
                solut_partial.cost = solut_crnt.cost
                iterILS = 0

            #t_perturb -= time.time()
            solut_crnt.s = perturb(solut_partial.s, info)
            #t_perturb += time.time()
            subseq_load(solut_crnt, info)
            iterILS += 1

        subseq_load(solut_partial, info)

        if solut_partial.cost < solut_best.cost:
            solut_best.s = solut_partial.s.copy()
            solut_best.cost = solut_partial.cost

        print("\tCurrent best solution cost: {}".format(solut_best.cost))

    print("COST: {}".format (solut_best.cost))
    print("SOLUTION: {}".format(solut_best.s))
   #print("Total Iterations RVND {}".format(IT))
   #print("swap", improv_swap)
   #print("two_opt", improv_two_opt)
   #print("reinsertion", improv_reinsertion)
   #print("or-opt_2", improv_or_2)
   #print("or-opt_3", improv_or_3)

def main() -> NoReturn:
    dimension :int
    cost : Matrix

    dimension , cost, rnd  = get_instance_info()

    R : List[float] = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]

    Imax : int = 10
    Iils : int = min(dimension, 100)

    info = tInfo(dimension, cost, rnd)

    """
    test = deepcopy(info)

    print(test.dimen)

    solut = tSolution([i for i in range(info.dimen)], subseq_fill(info), 0)

    solut.s.append(0)
    subseq_load(solut, info)

    print(solut.s, solut.seq, solut.cost)

"""

    start = time.time()
    GILS_RVND(Imax, Iils, R, info)
    print("TIME: %s " % (time.time() - start))

    print("ITERACOES: ", info.IT)

main()
"""
print("SWAP %s" % t_swap)
print("Reinsert %s" % t_reinsertion)
print("or_opt2 %s" % t_or_opt2)
print("or_opt3 %s" % t_or_opt3)
print("two_opt %s" % t_two_opt)
print("subseq_load %s" % t_seq)
"""
