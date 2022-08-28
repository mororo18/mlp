#! /usr/bin/bash
COMP=gfortran
rm *.o
$COMP -cpp -c data.f90
$COMP -cpp -c types.f90
$COMP -cpp -c assert.f90
$COMP -cpp -DOPa  main.f90 data.o types.o assert.f90 -o solve -g
