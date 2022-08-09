#! /usr/bin/bash

export DOTNET_PerfMapEnabled=1
export DOTNET_EnableDiagnostics=1
export DOTNET_EnableEventLog=1
export CORECLR_ENABLE_PROFILING=1

./build.sh
sudo python3 profile.py
