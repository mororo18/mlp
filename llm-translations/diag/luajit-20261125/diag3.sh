#!/usr/bin/env bash
# LuaJIT: tradução com a matriz de distâncias preenchida em ordem, sem buracos (como a nossa).
cd "$(dirname "$0")"; D=$PWD
P="ours_lua.literal_offsets cpp_lua.luajit_subseq_array cpp_lua.luajit_subseq_array_m tr_all_hoist_carry tr_all_hoist_carry_m"
for inst in gr202 kroA100; do
  cp $D/dm_$inst $D/w/distance_matrix
  for r in 1 2; do
    for d in $P; do
      out=$(cd $D/w/$d && luajit main.lua)
      t=$(echo "$out" | grep -E 'wall clock' | awk '{print $NF}'); [ -z "$t" ] && t=$(echo "$out" | grep -E '^TIME' | awk '{print $2}')
      echo "T $inst $d $r $t $(echo "$out" | grep COST | awk '{print $2}')"
    done
  done
done
echo FIM
