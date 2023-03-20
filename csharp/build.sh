#! /usr/bin/bash
# using mono
mcs -debug -optimize+ -out:solve_mcs *.cs 

# using dotnet
dotnet build
