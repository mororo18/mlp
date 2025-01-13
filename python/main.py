#! /usr/bin/pypy3

import time
from read import *
import math
import sys
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

    for i in range(info.dimen+1):
        arr.append([])
        for j in range(info.dimen+1):
            arr[i].append([0.0, 0.0, 0.0])

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

        i = int(len(c_list)*alpha) + 1

        sort(c_list, r, info)

        r_value = info.rnd[info.rnd_index]
        info.rnd_index += 1

        c = c_list[r_value]

        s.append(c)
        r = c
        c_list.remove(c)

    s.append(0)

    return s

def subseq_load(solut : tSolution, info : tInfo) -> NoReturn:
    for i in range(0, info.dimen+1):
        k : int = 1 - i - int(not i)

        solut.seq[i][i][tInfo.T] = 0.0
        solut.seq[i][i][tInfo.C] = 0.0
        solut.seq[i][i][tInfo.W] = float(int(not (i == 0)))

        for j in range(i+1, info.dimen+1):
        #while j < d:
            #print(i, j)
            j_prev : int = j - 1

            solut.seq[i][j][tInfo.T] = info.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][tInfo.T]
            solut.seq[i][j][tInfo.C] = solut.seq[i][j][tInfo.T] + solut.seq[i][j_prev][tInfo.C]
            solut.seq[i][j][tInfo.W] = float(j + k)

        #exit(0)

    solut.cost = solut.seq[0][info.dimen][tInfo.C] - tInfo.EPSILON

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

        cost_concat_1 = solut.seq[0][i_prev][tInfo.T] + info.cost[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[i][i_next][tInfo.T] + info.cost[solut.s[i]][solut.s[i_next+1]]

        cost_new = solut.seq[0][i_prev][tInfo.C]                                            + \
                solut.seq[i][i_next][tInfo.W]   * (cost_concat_1) + info.cost[solut.s[i_next]][solut.s[i]]  + \
                solut.seq[i_next+1][info.dimen][tInfo.W] * (cost_concat_2) + solut.seq[i_next+1][info.dimen][tInfo.C]

        if cost_new < cost_best:
            cost_best = cost_new - tInfo.EPSILON
            I = i
            J = i_next
            #print(cost_best, I, J)

        for j in range(i_next+1, (info.dimen)):
            j_next : int = j + 1
            j_prev : int = j - 1

            cost_concat_1 = solut.seq[0][i_prev][tInfo.T] + info.cost[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + info.cost[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][j_prev][tInfo.T] + info.cost[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3 + info.cost[solut.s[i]][solut.s[j_next]]

            cost_new = solut.seq[0][i_prev][tInfo.C]                                                + \
                    cost_concat_1                                                   + \
                    solut.seq[i_next][j_prev][tInfo.W] * cost_concat_2 + solut.seq[i_next][j_prev][tInfo.C] + \
                    cost_concat_3                                                   + \
                    solut.seq[j_next][info.dimen][tInfo.W] * cost_concat_4 + solut.seq[j_next][info.dimen][tInfo.C]

            if cost_new < cost_best:
                cost_best = cost_new - tInfo.EPSILON
                I = i
                J = j
                #print(cost_best, I, J)

    #print(cost_best, solut.cost, I, J)
    if cost_best < solut.cost - tInfo.EPSILON:
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
        
        reverse_cost : float = solut.seq[i][i+1][tInfo.T]
        for j in range(i+2, info.dimen):
            j_next : int = j + 1

            reverse_cost += info.cost[solut.s[j-1]][solut.s[j]] * (solut.seq[i][j][tInfo.W]-1)

            cost_concat_1 = solut.seq[0][i_prev][tInfo.T] + info.cost[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][j][tInfo.T] + info.cost[solut.s[j_next]][solut.s[i]]

            cost_new = solut.seq[0][i_prev][tInfo.C]                                        + \
                    solut.seq[i][j][tInfo.W] * cost_concat_1 + reverse_cost             + \
                    solut.seq[j_next][info.dimen][tInfo.W] * (cost_concat_2) + solut.seq[j_next][info.dimen][tInfo.C]

            #print(cost, i, "        ", j)
            #print(seq[C][0][i_prev], seq[W][i][j] * cost_concat_1,seq[W][j_next][n] * (cost_concat_2), seq[C][j_next][n] )
            if cost_new < cost_best:
                cost_best = cost_new - tInfo.EPSILON
                I = i
                J = j

            #if i == 2 and j == 10:
                #print ("opaa")

    if cost_best < solut.cost - tInfo.EPSILON:
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


            cost_concat_1 =                 seq[0]     [k]     [tInfo.T] + info.cost[solut.s[k]]     [solut.s[i]]
            cost_concat_2 = cost_concat_1 + seq[i]     [j]     [tInfo.T] + info.cost[solut.s[j]]     [solut.s[k_next]]
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev][tInfo.T] + info.cost[solut.s[i_prev]][solut.s[j_next]]


            cost_new = seq[0][k][tInfo.C]                                                                     + \
                    seq[i]     [j]         [tInfo.W]  * cost_concat_1 + seq[i]     [j]         [tInfo.C] + \
                    seq[k_next][i_prev]    [tInfo.W]  * cost_concat_2 + seq[k_next][i_prev]    [tInfo.C] + \
                    seq[j_next][info.dimen][tInfo.W]  * cost_concat_3 + seq[j_next][info.dimen][tInfo.C]

            if cost_new < cost_best:
                cost_best = cost_new - tInfo.EPSILON
                I = i
                J = j
                POS = k


        for k in range(i+opt, info.dimen):
            k_next = k+1

            cost_concat_1 =                 seq[0]     [i_prev][tInfo.T] + info.cost[solut.s[i_prev]][solut.s[j_next]]
            cost_concat_2 = cost_concat_1 + seq[j_next][k]     [tInfo.T] + info.cost[solut.s[k]]     [solut.s[i]]
            cost_concat_3 = cost_concat_2 + seq[i]     [j]     [tInfo.T] + info.cost[solut.s[j]]     [solut.s[k_next]]


            cost_new = seq[0][i_prev][tInfo.C]                                        + \
                    seq[j_next][k][tInfo.W]   * cost_concat_1 + seq[j_next][k][tInfo.C] + \
                    seq[i][j][tInfo.W]        * cost_concat_2 + seq[i][j][tInfo.C]      + \
                    seq[k_next][info.dimen][tInfo.W]   * cost_concat_3 + seq[k_next][info.dimen][tInfo.C]

            if cost_new < cost_best:
                cost_best = cost_new - tInfo.EPSILON
                I = i
                J = j
                POS = k

    if cost_best < solut.cost - tInfo.EPSILON:
        reinsert(solut.s, I, J, POS+1)
        subseq_load(solut, info)
        """
        print("\nreinsertion", opt)
        print(cost_best)
        print(solut.cost)
        """
        return True

    return False

def RVND(solut : tSolution, info : tInfo) -> NoReturn:

    neighbd_list : List[int] = [tInfo.SWAP, tInfo.TWO_OPT, tInfo.REINSERTION, tInfo.OR_OPT_2, tInfo.OR_OPT_3]

    while len(neighbd_list) > 0:
        i = info.rnd[info.rnd_index]
        info.rnd_index += 1


        neighbd = neighbd_list[i]

        improve : bool = False

        if neighbd == tInfo.SWAP:
            improve = search_swap(solut, info)
        elif neighbd == tInfo.TWO_OPT:
            improve = search_two_opt(solut, info)
        elif neighbd == tInfo.REINSERTION:
            improve = search_reinsertion(solut, info, tInfo.REINSERTION, solut.seq)
        elif neighbd == tInfo.OR_OPT_2:
            improve = search_reinsertion(solut, info, tInfo.OR_OPT_2, solut.seq)
        elif neighbd == tInfo.OR_OPT_3:
            improve = search_reinsertion(solut, info, tInfo.OR_OPT_3, solut.seq)

        if improve == True:
            neighbd_list = [tInfo.SWAP, tInfo.TWO_OPT, tInfo.REINSERTION, tInfo.OR_OPT_2, tInfo.OR_OPT_3]
            
        else:
            neighbd_list.pop(i)

def perturb(sl : List[int], info : tInfo) -> List[int]:
    s = sl.copy()

    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

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

def main() -> NoReturn:
    dimension :int
    cost : Matrix

    dimension , cost, rnd  = get_instance_info()

    R : List[float] = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]

    Imax : int = 10
    Iils : int = min(dimension, 100)

    info = tInfo(dimension, cost, rnd)


    start = time.time()
    GILS_RVND(Imax, Iils, R, info)
    print("TIME: %s " % (time.time() - start))

    print("ITERACOES: ", info.IT)

main()
