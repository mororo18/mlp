#! /usr/bin/bash
g++ -O3 -ffast-math -o mlp main.cpp Data.cpp && ./mlp "$@"
