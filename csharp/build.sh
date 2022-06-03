#! /bin/bash
# using mono
mcs -optimize+ main.cs Data.cs GILS_RVND.cs  

# using dotnet
dotnet build --output ./dotnet_output
