#! /usr/bin/bash
# using mono
mcs -optimize+ -out:solve *.cs 

# using dotnet
rm bin/Debug/net6.0/csharp
dotnet build
