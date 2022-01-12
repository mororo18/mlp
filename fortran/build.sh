#! /usr/bin/bash
gfortran -c Data.f90
gfortran -O3 main.f90 Data.o -o main -g 
