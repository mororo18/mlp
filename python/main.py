#! /usr/bin/python3

from test import *
from random import randint

n, m = get_instance_info()

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

    print(s)

    return s

def RVND(s):
    return None

def perturb(s):
    return None

def cost_calc(s):
    cost = 0
    # calc errado
    for i in range(len(s) - 1):
        cost += m[s[i]][s[i+1]]
        
    return cost

def GILS_RVND(Imax, Iils, R):

    cost_best = float('inf')

    for i in range(Imax):
        alpha = R[randint(0, len(R)-1)]

        s = construction(alpha)
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
