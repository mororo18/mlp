#!/usr/bin/env bash
# C# na berlin52: qual parte do JIT em camadas pesa, com e sem a busca de aquecimento.
cd "$(dirname "$0")"; D=$PWD
run(){ local d=$1 dll=$2; shift 2
  out=$(cd $D/w/$d && env "$@" dotnet $dll)
  echo "$out" | grep -E 'wall clock' | awk '{print $NF}' | tr , .; }
cp $D/dm_berlin52 $D/w/distance_matrix
for r in 1 2 3 4; do
  for cfg in "default:X=1" "tc0:DOTNET_TieredCompilation=0" "qjfl0:DOTNET_TC_QuickJitForLoops=0" "osr0:DOTNET_TC_OnStackReplacement=0" "ccd0:DOTNET_TC_CallCountingDelayMs=0" "bg0:DOTNET_TC_BackgroundWorkerTimeoutMs=0" "qj0:DOTNET_TC_QuickJit=0"; do
    n=${cfg%%:*}; e=${cfg#*:}
    echo "T berlin52 ours $n $r $(run cs_ours bin/Release/net8.0/csharp.dll $e)"
    echo "T berlin52 tr $n $r $(run cs_tr bin/Release/net8.0/solve.dll $e)"
    if [ $n = default ] || [ $n = tc0 ] || [ $n = qjfl0 ]; then
      echo "T berlin52 ours_warm $n $r $(run cs_ours_warm bin/Release/net8.0/csharp.dll $e)"
      echo "T berlin52 tr_warm $n $r $(run cs_tr_warm bin/Release/net8.0/solve.dll $e)"
    fi
  done
done
cp $D/dm_kroA100 $D/w/distance_matrix
for r in 1 2; do
  for cfg in "tc0:DOTNET_TieredCompilation=0" "qjfl0:DOTNET_TC_QuickJitForLoops=0"; do
    n=${cfg%%:*}; e=${cfg#*:}
    echo "T kroA100 ours $n $r $(run cs_ours bin/Release/net8.0/csharp.dll $e)"
    echo "T kroA100 tr $n $r $(run cs_tr bin/Release/net8.0/solve.dll $e)"
  done
done
echo FIM
