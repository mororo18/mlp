#! /bin/bash
rm main
gfortran -c Data.f90
gfortran -O3 main.f90 Data.o -o main -g 
