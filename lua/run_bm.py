import os

for i in range(10):
    os.system("luajit -O8 main.lua")
