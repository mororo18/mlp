#! /usr/bin/env python3


def matrix_fill(n):
    matrix = []
    for i in range(n):
        line = []
        for j in range(n):
            line.append(float('inf'))

        matrix.append(line)


    return matrix

def get_instance_info(f_name='../distance_matrix'):
    source_f = open(f_name, 'r')
    cost = []

    i = -1  
    dimension = None

    for line in source_f:
        if i == -1:
            dimension = int(line)
            cost = matrix_fill(dimension)

        else:
            j = i+1
            while(line.find(' ') != -1):
                index = line.find(' ')

                cost[i][j] = float(line[:index])
                cost[j][i] = cost[i][j]

                j += 1
                line = line[index+1:]

            #print(cost[i])

        i += 1

    #print(dimension)
    #exit(1)

    return dimension, cost
