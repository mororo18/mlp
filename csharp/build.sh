#! /usr/bin/bash
# using mono
mcs -optimize+ -out:solve *.cs 

# using dotnet
dotnet build
