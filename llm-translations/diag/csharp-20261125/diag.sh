#!/usr/bin/env bash
# Roda dentro do container nix (dotnet da rodada), preso a um núcleo.
# C# na berlin52: nossa e tradução (S1+S2) sob variáveis do runtime (GC e JIT em camadas).
cd "$(dirname "$0")"; D=$PWD
dotnet --version
run(){ # dir dll envs...
  local d=$1 dll=$2; shift 2
  out=$(cd $D/w/$d && env "$@" dotnet $dll)
  echo "$out" | grep -E 'wall clock' | awk '{print $NF}' | tr , .
}
for inst in berlin52 kroA100; do
  cp $D/dm_$inst $D/w/distance_matrix
  R=5; [ $inst = kroA100 ] && R=2
  for r in $(seq 1 $R); do
    for cfg in "default:X=1" "gcServer1:DOTNET_gcServer=1" "gcServer0:DOTNET_gcServer=0" "gcConc0:DOTNET_gcConcurrent=0" "pgo0:DOTNET_TieredPGO=0" "tc0:DOTNET_TieredCompilation=0"; do
      n=${cfg%%:*}; e=${cfg#*:}
      [ $inst = kroA100 ] && [ $n != default ] && [ $n != gcServer1 ] && [ $n != gcServer0 ] && continue
      echo "T $inst ours $n $r $(run cs_ours bin/Release/net8.0/csharp.dll $e)"
      echo "T $inst tr $n $r $(run cs_tr bin/Release/net8.0/solve.dll $e)"
    done
  done
done
echo FIM
