#! /usr/bin/python3

from read import *
from random import randint
#from rcdtype import *
import types

n, m = get_instance_info()

#subseq_info = types.SimpleNamespace('T C W')

def subseq_info_fill(n):
    matrix = []

    for i in range(n+1):
        matrix.append([])
        for j in range(n+1):
            matrix[i].append(types.SimpleNamespace())
            matrix[i][j].T = 0
            matrix[i][j].C = 0
            matrix[i][j].W = 0
            #matrix[i].append(subseq_info(0, 0, 0))

        print(matrix[i])



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

        seq[i][i].T = 0
        seq[i][i].C = 0
        seq[i][i].W = int(not (i == 0))

        j = i + 1
        while j < d:
            a = j - 1

            seq[i][j].T = m[s[a]][s[j]] + seq[i][a].T
            seq[i][j].C = seq[i][j].T + seq[i][a].C
            seq[i][j].W = j + k

            j += 1

        #exit(1)

        i += 1

def RVND(s):
    return None

def perturb(s):
    return None

def GILS_RVND(Imax, Iils, R):

    cost_best = float('inf')


    subseq = subseq_info_fill(n)

    for i in range(Imax):
        alpha = R[randint(0, len(R)-1)]

        s = construction(alpha)
        #s = [0, 7, 8, 10, 12, 6, 11, 5, 4, 3, 2, 13, 1, 9, 0]
        # s_cost igual a 20315 p burma14
        print('aqui')
        subseq_info_load(s, subseq)
        print(subseq[0][n].C)
        exit(1)
        sl = s

        iterILS = 0
        while iterILS < Iils:
            s = RVND(s)

            if cost_calc(s) < cost_calc(sl):
                sl = s
                iterILS = 0

            s = perturb(sl)
            iterILS += 1


def main():
    
    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]
    
    Imax = 10
    Iils = min(n, 100)

    GILS_RVND(Imax, Iils, R)


main()
