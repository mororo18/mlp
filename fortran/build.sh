#! /bin/bash
rm main
gfortran -c Data.f90
gfortran -Ofast main.f90 Data.o -o main -g 
