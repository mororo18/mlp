#!/usr/bin/env bash
# LuaJIT do ambiente nix: tradução (L1, L1-L4, L1-L5) e nossa com LJ1, gr202 e kroA100; -jv no gr202.
cd "$(dirname "$0")"; D=$PWD
luajit -v; which luajit
P="ours_lua.literal_offsets cpp_lua.luajit_subseq_array tr_all_hoist tr_all_hoist_carry"
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
cp $D/dm_gr202 $D/w/distance_matrix
for d in ours_lua.literal_offsets cpp_lua.luajit_subseq_array; do
  (cd $D/w/$d && luajit -jv=$D/jv_gr202_$d.txt main.lua > /dev/null)
done
echo FIM
