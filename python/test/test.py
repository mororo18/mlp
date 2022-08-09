import timeit

A = open("A.py").read()
#B = open("B.py").read()
#C = open("C.py").read()

ar = timeit.timeit("search(np_arr, n)", setup=A)
dt = timeit.timeit("search(dicti, n)", setup=A)
lt = timeit.timeit("search(lista, n)", setup=A)

print("np_array", ar)
print("dict", dt)
print("list", lt)

#print(" {} times".format(ar / dt))
