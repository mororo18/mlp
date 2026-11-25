#!/usr/bin/env bash
# Roda dentro do container nix (lua 5.4 da rodada), preso a um núcleo.
cd "$(dirname "$0")"; D=$PWD
lua -v
for inst in berlin52 kroA100; do
  cp $D/dm_$inst $D/w/distance_matrix
  R=5; [ $inst = kroA100 ] && R=2
  for r in $(seq 1 $R); do
    for d in lua_ours lua_z lua_tr; do
      out=$(cd $D/w/$d && lua main.lua)
      t=$(echo "$out" | grep -E 'wall clock' | awk '{print $NF}'); [ -z "$t" ] && t=$(echo "$out" | grep -E '^TIME' | awk '{print $2}')
      echo "T $inst $d $r $t $(echo "$out" | grep COST | awk '{print $2}')"
    done
  done
done
echo FIM
