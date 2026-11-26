#!/usr/bin/env bash
# Transplante (swap, 2-opt, os dois), D7 (resultado por referência) e D8 (matriz em struct),
# binários estáticos do GCC 14 no host, núcleo $1. Espera a diag-st2 terminar.
cd "$(dirname "$0")"; C=${1:-16}
until grep -q "^FIM" ../diag-st2/out.txt 2>/dev/null; do sleep 15; done
P="ours tr t_swap t_2opt t_both d7_ref d8_struct"
for inst in gr202 lin318; do
  cp dm_$inst w/distance_matrix
  for r in 1 2 3; do for v in $P; do
    out=$(cd w/$v && taskset -c $C ./solve)
    echo "T $inst $v $r $(echo "$out" | grep -E 'wall clock' | awk '{print $NF}') $(echo "$out" | grep -o 'COST: *[0-9]*')"
  done; done
done
cat /proc/loadavg
echo FIM
