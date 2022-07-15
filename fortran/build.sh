#! /usr/bin/bash
COMP=gfortran
$COMP -c Data.f90
$COMP -Ofast main.f90 Data.o -o solve -g
