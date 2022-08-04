#! /usr/bin/bash
# using mono
mcs -optimize+ -out:solve main.cs Data.cs GILS_RVND.cs  

# using dotnet
dotnet build
