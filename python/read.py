#! /usr/bin/env python3


def matrix_fill(n:int):
    matrix = []
    for i in range(n):
        line = []
        for j in range(n):
            line.append(float('inf'))

        matrix.append(line)


    return matrix

def get_instance_info(f_name='../distance_matrix'):
    f = open(f_name, 'r')
    cost = []

    i = -1  

    dimension = int(f.readline())
    cost = matrix_fill(dimension)

    for i in range(dimension):
        line = f.readline()
        j = i+1
        while(line.find(' ') != -1):
            index = line.find(' ')

            cost[i][j] = float(line[:index])
            cost[j][i] = cost[i][j]

            j += 1
            line = line[index+1:]

    f.readline()
    f.readline()
    rnd_size = int(f.readline())
    rnd = []
    for i in range(rnd_size):
        rnd_value = int(f.readline())
        rnd.append(rnd_value)

    return dimension, cost, rnd
