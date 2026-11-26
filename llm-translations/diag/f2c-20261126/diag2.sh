#!/usr/bin/env bash
# Ciclos (rdtsc) por busca, fora dos laços: reinserção, Or-opt2, Or-opt3, swap, 2-opt e atualização.
cd "$(dirname "$0")"; D=$PWD
F="-m64 -Ofast -fPIC -fexceptions -DNDEBUG -DIL_STD -Winline -g -std=c++0x"
for v in ours_t int_t; do (cd w/$v && g++ $F -c main.cpp -o main.o 2>/dev/null && g++ $F -c readData.cpp -o readData.o 2>/dev/null && g++ $F main.o readData.o -o solve) || echo FAIL $v; done
(cd w/tr_t && g++ $F -DDOUBLE -c main.cpp -o main.o 2>/dev/null && g++ $F -DDOUBLE main.o -o solve) || echo FAIL tr_t
for inst in gr202 lin318; do
  cp $D/dm_$inst $D/w/distance_matrix
  for r in 1 2; do
    for v in ours_t int_t tr_t; do
      out=$(cd w/$v && ./solve 2>&1)
      echo "C $inst $v $r $(echo "$out" | grep -E 'wall clock' | awk '{print $NF}') $(echo "$out" | grep '^CYC' | cut -d' ' -f2-)"
    done
  done
done
echo FIM
