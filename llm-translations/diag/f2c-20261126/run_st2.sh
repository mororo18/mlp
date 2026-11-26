#!/usr/bin/env bash
# Binários estáticos (GCC 14, notebook), um por diferença, rodando direto no host, núcleo $1.
cd "$(dirname "$0")"; C=${1:-16}
P="ours int i64 tr d2_vec d3_1tu d4_noalign d5_ours_al d5_tr_al"
for inst in gr202 lin318; do
  cp dm_$inst w/distance_matrix
  for r in 1 2 3; do for v in $P; do
    out=$(cd w/$v && taskset -c $C ./solve)
    echo "T $inst $v $r $(echo "$out" | grep -E 'wall clock' | awk '{print $NF}') $(echo "$out" | grep -o 'COST: *[0-9]*')"
  done; done
done
cat /proc/loadavg
echo FIM
