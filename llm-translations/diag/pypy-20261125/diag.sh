#!/usr/bin/env bash
# Roda dentro do container nix (pypy3 da rodada), preso a um núcleo.
cd "$(dirname "$0")"; D=$PWD; CPU=${CPU:-14}
pypy3 --version 2>&1 | tail -1
prog(){ case $1 in py_tr*) echo mlp_solver.py;; *) echo main_pypy.py;; esac; }
for inst in berlin52 kroA100; do
  cp $D/dm_$inst $D/w/distance_matrix
  for r in 1 2 3 4 5 6; do
    for d in py_ours py_tr py_ours_warm py_tr_warm; do
      out=$(cd $D/w/$d && pypy3 $(prog $d))
      t=$(echo "$out" | grep -E 'wall clock' | awk '{print $NF}'); [ -z "$t" ] && t=$(echo "$out" | grep -E '^TIME' | awk '{print $2}')
      echo "T $inst $d $r $t $(echo "$out" | grep COST | awk '{print $2}')"
    done
  done
  for d in py_ours py_tr py_ours_warm py_tr_warm; do
    (cd $D/w/$d && PYPYLOG=jit-summary:$D/jit_${inst}_$d.log pypy3 $(prog $d) > /dev/null)
  done
done
echo FIM
