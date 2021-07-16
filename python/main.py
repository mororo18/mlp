#! /usr/bin/python3

from read import *
from random import randint
#from rcdtype import *
import types

EPSILON = 0.000000001

SWAP        = 0
TWO_OPT     = 1
REINSERTION = 2
OR_OPT_2    = 3
OR_OPT_3    = 4

n, m = get_instance_info()

W = "W"
T = "T"
C = "C"

improv_flag = False

#subseq_info = types.SimpleNamespace('T C W')

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
            #matrix[i].append(types.SimpleNamespace())
            matrix[W][i].append(0.0)
            matrix[T][i].append(0.0)
            matrix[C][i].append(0.0)
            #matrix[i].append(subseq_info(0, 0, 0))



    '''
    print(matrix[W])
    print(matrix[T])
    print(matrix[C])
    exit(1)
    '''

    return matrix

def construction(alpha):
    s = [0]
    c_list = []
    for i in range(1, n):
        c_list.append(i)

    r = 0

    while len(c_list) > 0:
        c_list = sorted(c_list, key = lambda i : m[i][r], reverse=False)
        print(c_list)

        '''
        for i in range(len(c_list)):
            print(m[c_list[i]][r])
        '''

        i = int(len(c_list)*alpha)
        if i == 0:
            Rc_list = [c_list[i]]
        else:
            Rc_list = c_list[:i]

        print(Rc_list)

        c = Rc_list[randint(0, len(Rc_list)-1)]
        s.append(c)
        r = c
        c_list.remove(r)

    s.append(0)

    print(s)

    return s

def subseq_info_load(s, seq):
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

            '''
            print(s[a], s[j], a, j)
            print(type(m[s[a]][s[j]]), type(seq[T][i][a]))
            print(m[s[a]][s[j]],  seq[T][i][a])
            '''
            seq[T][i][j] = m[s[a]][s[j]] + seq[T][i][a]
            seq[C][i][j] = seq[T][i][j] + seq[C][i][a]
            seq[W][i][j] = j + k

            seq[T][j][i] = seq[T][i][j]
            seq[C][j][i] = seq[C][i][j]
            seq[W][j][i] = seq[W][i][j]

            j += 1

        #print(seq[C][i])

        i += 1
    #exit(1)

def swap(s, i, j):
    s[i], s[j] = s[j], s[i]

def reverse(s, i, j):
    s[i:j+1] = s[i:j+1][::-1]

def reinsert(s, i, j, pos):
    if i < pos:
        s[pos:pos] = s[i:j+1]
        s[:-1] = s[:i] + s[j+1:]
    else:
        sub = s[i:j+1]
        s[:-1] = s[:i] + s[j+1:]
        s[pos:pos] = sub


def search_swap(s, seq):
    cost_best = float('inf')
    I = None
    J = None

    for i in range(1, (n-1) ):
        i_prev = i - 1
        i_next = i + 1

        total = seq[T][0][i_prev] + m[s[i_prev]][s[i_next]]

        cost = seq[C][0][i_prev] + seq[W][i][i_next] * (total) + m[s[i_next]][s[i]]
        cost += seq[W][i_next+1][n] * (total + seq[T][i][i_next] + m[s[i]][s[i_next+1]]) + seq[C][i_next+1][n]

        if cost < cost_best:
            cost_best = cost - EPSILON
            I = i
            J = i_next

        #if i == n-2:
        #    continue

        for j in range(i_next+1, (n)):
            j_next = j + 1
            j_prev = j - 1

            #print(seq[T][0][i_prev])
            #print( m[s[i_prev]][s[j]], s[i_prev], s[j], i_prev, j,  "\n")
            total_1 = seq[T][0][i_prev] + m[s[i_prev]][s[j]]
            #TODO problema aqui
            print(total_1, m[s[j]][s[i_next]], s[j], s[i_next], j, i_next)
            total_2 = total_1 + m[s[j]][s[i_next]]
            total_3 = total_2 + seq[T][i_next][j_prev] + m[s[j_prev]][s[i]]
            #TODO problema aqui
            print(total_3, m[s[i]][s[j_next]], s[i], s[j_next], i, j_next)
            total_4 = total_3 + m[s[i]][s[j_next]]

            cost = seq[C][0][i_prev] + total_1 + seq[W][i_next][j_prev] * total_2 + seq[C][i_next][j_prev] + total_3 + seq[W][j_next][n] * total_4 + seq[C][j_next][n]
            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j


    if cost_best < seq[C][0][n] - EPSILON:
        swap(s, I, J)
        print("swap", I, J)
        print(s)
        subseq_info_load(s, seq)
        improv_flag = True

    #swap(s, 7, 4)
    #return None

def search_two_opt(s, seq):
    cost_best = float('inf')

    for i in range(1, n-1):
        i_prev = i - 1
        for j in range(i+2, n):
            j_next = j + 1

            total = seq[T][0][i_prev] + m[s[j]][s[i_prev]]

            cost = seq[C][0][i_prev] + seq[W][i][j] * total + sum([seq[T][x][j] for x in range(i, j)])
            #TODO problema aqui
            cost += seq[W][j_next][n] * (total + seq[T][i][j] + m[s[j_next]][s[i]]) + seq[C][j_next][n]

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j

    if cost_best < seq[C][0][n] - EPSILON:
        reverse(s, I, J)
        print("reverse", I, J)
        print(s)
        subseq_info_load(s, seq)
        improv_flag = True

def search_reinsertion(s, seq, OPT):
    cost_best = float('inf')
    I = None
    J = None
    POS = None
    opt = OPT - 1 

    for i in range(1, n - opt + 1):
        j = opt + i - 1

        j_next = j+1
        i_prev = i-1

        for k in range(0, i_prev):
            k_next = k+1

            total_1 = seq[T][0][k] + m[s[k]][s[i]]
            total_2 = total_1 + seq[T][i][j] + m[s[j]][s[k_next]]
            
            cost = seq[C][0][k] + seq[W][i][j] * total_1 + seq[C][i][j]
            cost += seq[W][k_next][i_prev] * total_2 + seq[C][k_next][i_prev]
            cost += seq[W][j_next][n] * (total_2 + seq[T][k_next][i_prev] + m[s[i_prev]][s[j_next]]) + seq[C][j_next][n]

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j
                POS = k

        for k in range(i+opt, n - opt - 1):
            k_next = k+1
            
            total_1 = seq[T][0][i_prev] + m[s[i_prev]][s[j_next]]
            #TODO problema aqui
            total_2 = total_1 + seq[T][j_next][k] + m[s[k]][s[i]]
            
            cost = seq[C][0][i_prev] + seq[W][j_next][k] * total_1 + seq[C][j_next][k]
            cost += seq[W][i][j] * total_2 + seq[C][i][j]
            #TODO problema aqui
            cost += seq[W][k_next][n] * (total_2 + seq[T][i][j] + m[s[j]][s[k_next]]) + seq[C][k_next][n]

            if cost < cost_best:
                cost_best = cost - EPSILON
                I = i
                J = j
                POS = k



    if cost_best < seq[C][0][n] - EPSILON:
        reinsert(s, I, J, POS+1)
        subseq_info_load(s, seq)
        improv_flag = True
        print("reinsertion", I,"a", J, "em" , POS)
        print(s)
    #return None

def RVND(s, subseq):

    #neighbd_list = [SWAP, TWO_OPT]
    #neighbd_list = [TWO_OPT]
    #neighbd_list = [REINSERTION]
    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]

    while len(neighbd_list) > 0:
        i = randint(0, len(neighbd_list)-1)
        neighbd = neighbd_list[i]

        improv_flag = False

        if neighbd == SWAP:
            search_swap(s, subseq)
        elif neighbd == TWO_OPT:
            search_two_opt(s, subseq)
        elif neighbd == REINSERTION:
            search_reinsertion(s, subseq, REINSERTION)
        elif neighbd == OR_OPT_2:
            search_reinsertion(s, subseq, OR_OPT_2)
        elif neighbd == OR_OPT_3:
            search_reinsertion(s, subseq, OR_OPT_3)


        if improv_flag == True:
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]
        else:
            neighbd_list.pop(i)
        

    print(subseq[C][0][n])
    #exit(1)

    return s

def perturb(s):
    return None

def GILS_RVND(Imax, Iils, R):

    cost_best = float('inf')
    s_best = []

    subseq = subseq_info_fill(n)

    for i in range(Imax):
        alpha = R[randint(0, len(R)-1)]

        s = construction(alpha)
        #s = [0, 7, 8, 10, 12, 6, 11, 5, 4, 3, 2, 13, 1, 9, 0]
        # s_cost igual a 20315 p burma14
        print('aqui')
        subseq_info_load(s, subseq)
        print(subseq[C][0][n])
        #exit(1)
        sl = s
        rvnd_cost_best = subseq[C][0][n] - EPSILON

        iterILS = 0
        while iterILS < Iils:
            s = RVND(s, subseq)
            rvnd_cost_crnt  = subseq[C][0][n] - EPSILON
            if rvnd_cost_crnt < rvnd_cost_best:
                rvnd_cost_best = rvnd_cost_crnt
                sl = s
                iterILS = 0

            #s = perturb(sl)
            subseq_info_load(s, subseq)
            iterILS += 1

        subseq_info_load(sl, subseq)
        sl_cost = subseq[C][0][n] - EPSILON

        if sl_cost < cost_best:
            s_best = sl
            cost_best = sl_cost

    print(cost_best)




def main():
    
    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]
    
    Imax = 10
    Iils = min(n, 100)

    GILS_RVND(Imax, Iils, R)


main()
