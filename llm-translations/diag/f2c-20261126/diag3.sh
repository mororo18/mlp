#!/usr/bin/env bash
# cachegrind (GCC 15, binários de diag.sh) no kroA100: instruções e cache simulada por função.
cd "$(dirname "$0")"; D=$PWD
mkdir -p $D/cg; cp $D/dm_kroA100 $D/cg/distance_matrix
for v in ours int tr; do
  mkdir -p $D/cg/$v; cp $D/w/$v/solve $D/cg/$v/
  (cd $D/cg/$v && valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=$D/cg/$v.out ./solve > $D/cg/$v.stdout 2> $D/cg/$v.stderr)
  cg_annotate --show=Ir,D1mr,DLmr $D/cg/$v.out 2>/dev/null | head -60 > $D/cg/$v.txt
done
echo FIM
