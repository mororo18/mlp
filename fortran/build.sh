#! /bin/bash
gfortran-7 -c Data.f90
gfortran-7 -O3 main.f90 Data.o -o main -g 
