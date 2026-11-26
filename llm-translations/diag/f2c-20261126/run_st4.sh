#!/usr/bin/env bash
# Combinações D8 + swap transplantado + índice int no nosso C++ (binários estáticos do GCC 14),
# no host, núcleo $1, 3 execuções alternadas em gr202 e lin318.
cd "$(dirname "$0")"; C=${1:-16}
P="ours tr int d8_struct t_swap c_d8_int c_d8_swap c_all"
for inst in gr202 lin318; do
  cp dm_$inst w/distance_matrix
  for r in 1 2 3; do for v in $P; do
    out=$(cd w/$v && taskset -c $C ./solve)
    echo "T $inst $v $r $(echo "$out" | grep -E 'wall clock' | awk '{print $NF}') $(echo "$out" | grep -o 'COST: *[0-9]*')"
  done; done
done
cat /proc/loadavg
echo FIM
