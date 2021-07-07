#! /usr/bin/python3

from test import *
from random import randint

n, m = get_instance_info()

def construct(alpha):
    return

def RVND(s):
    return

def perturb(s):
    return

def cost_calc(s):
    cost = 0
    for i in range(len(s) - 1):
        cost += m[s[i]][s[i+1]]
        
    return cost

def GILS_RVND(Imax, Iils, R):

    cost_best = float('inf')

    for i in range(Imax):
        alpha = R[randint(0, len(R))]

        s = construction(alpha)
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



