#! /usr/bin/bash
gfortran -c Data.f90
gfortran main.f90 Data.o -o main
