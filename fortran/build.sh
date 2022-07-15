#! /usr/bin/bash
gfortran -c Data.f90
gfortran -Ofast main.f90 Data.o -o solve -g
