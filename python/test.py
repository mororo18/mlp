import timeit

A = open("A.py").read()
B = open("B.py").read()
C = open("C.py").read()

ar = timeit.timeit("a[0]", setup=A)
dt = timeit.timeit("a[0]", setup=B)
lt = timeit.timeit("a[0]", setup=C)

print("array", ar)
print("dict", dt)
print("list", lt)

print(" {} times".format(ar / dt))
