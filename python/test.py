from struct import *

l = []

for i in range(1, 10):
    l.append(pack('id', i, 0.00034))
    print(pack('i', i))

for i in range(1, 10):
    print(unpack('id', l[i-1]))
