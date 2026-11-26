#!/usr/bin/env bash
# Diferenças D2 a D5 (GCC 15), cada uma sozinha no nosso C++, contra o nosso, o índice int e a tradução.
cd "$(dirname "$0")"; D=$PWD
F="-m64 -Ofast -fPIC -fexceptions -DNDEBUG -DIL_STD -Winline -g -std=c++0x"
AL="-falign-loops=32 -falign-jumps=32 -falign-functions=64"
b2(){ (cd w/$1 && g++ $F $2 -c main.cpp -o main.o 2>/dev/null && g++ $F $2 -c readData.cpp -o readData.o 2>/dev/null && g++ $F $2 main.o readData.o -o solve) || echo FAIL $1; }
b2 d2_vec ""; b2 d4_noalign ""; b2 d5_ours_al "$AL"
(cd w/d3_1tu && g++ $F -c main.cpp -o main.o 2>/dev/null && g++ $F main.o -o solve) || echo FAIL d3_1tu
(cd w/d5_tr_al && g++ $F $AL -DDOUBLE -c main.cpp -o main.o 2>/dev/null && g++ $F $AL -DDOUBLE main.o -o solve) || echo FAIL d5_tr_al
P="ours int tr d2_vec d3_1tu d4_noalign d5_ours_al d5_tr_al"
for inst in gr202 lin318; do
  cp $D/dm_$inst $D/w/distance_matrix
  for r in 1 2 3; do for v in $P; do
    out=$(cd w/$v && ./solve)
    echo "T $inst $v $r $(echo "$out" | grep -E 'wall clock' | awk '{print $NF}') $(echo "$out" | grep -o 'COST: *[0-9]*')"
  done; done
done
cat /proc/loadavg
echo FIM
