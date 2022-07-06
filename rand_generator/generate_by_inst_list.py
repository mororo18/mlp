import os
f = open("../inst-list")

os.chdir("../")

for inst in f:
    os.system("./load " + inst)
    os.chdir("rand_generator")
    os.system("julia main.jl")
    os.chdir("../")
    os.system("./load " + inst)
