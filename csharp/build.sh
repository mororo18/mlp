#! /usr/bin/bash
# using mono
mcs -debug -optimize+ -out:solve *.cs 

# using dotnet
dotnet build
