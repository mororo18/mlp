import numpy as np
import array 
from random import randint
n = 10

index = []
for i in range(n):
    index.append(randint(0, n-1))


# Init LIST
lista = []
for i in range(n):
    lista.append(list(range(n)))

# Init npArray
np_arr = np.zeros([n, n])


dicti = {}
for i in range(n):
    dicti[i] = dict.fromkeys(range(n))
    for j in range(n):
        dicti[i][j] = 0.0


#arrList = [array.array('f', list(range(n))) for i in range(n)]
##arr = array.array('d', arrList)

def search(arr, n):
    for i in range(n):
        for j in range(n):
            value = arr[index[i]][index[j]] + arr[i][j]

