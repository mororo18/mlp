#!/usr/bin/env bash
# Grade de parâmetros do JIT do PyPy (diagnóstico do aquecimento).
cd "$(dirname "$0")"; D=$PWD
prog(){ case $1 in py_tr*) echo mlp_solver.py;; *) echo main_pypy.py;; esac; }
CFGS="default threshold=200,function_threshold=300 threshold=50,function_threshold=80 trace_eagerness=50 threshold=200,function_threshold=300,trace_eagerness=50 trace_limit=20000"
for inst in berlin52 kroA100; do
  cp $D/dm_$inst $D/w/distance_matrix
  for r in 1 2 3; do
    for cfg in $CFGS; do
      J=""; [ "$cfg" != default ] && J="--jit $cfg"
      for d in py_ours py_ours_pp1 py_tr; do
        out=$(cd $D/w/$d && pypy3 $J $(prog $d))
        t=$(echo "$out" | grep -E 'wall clock' | awk '{print $NF}'); [ -z "$t" ] && t=$(echo "$out" | grep -E '^TIME' | awk '{print $2}')
        echo "T $inst $cfg $d $r $t $(echo "$out" | grep COST | awk '{print $2}')"
      done
    done
  done
done
echo FIM
