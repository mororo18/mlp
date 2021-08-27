# !/usr/bin/bash
dotnet build --output ./build_output
dotnet ./build_output/csharp.dll
