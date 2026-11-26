#!/usr/bin/env bash
# Fortran -> C++ no GCC 15 do ambiente Nix: build com as flags do harness, código de máquina das
# funções quentes e tempos alternados no gr202 e no lin318.
cd "$(dirname "$0")"; D=$PWD
g++ --version | head -1
F="-m64 -Ofast -fPIC -fexceptions -DNDEBUG -DIL_STD -Winline -g -std=c++0x"
for v in ours int i64; do (cd w/$v && g++ $F -c main.cpp -o main.o 2>/dev/null && g++ $F -c readData.cpp -o readData.o 2>/dev/null && g++ $F main.o readData.o -o solve && g++ $F -S -fverbose-asm main.cpp -o main.s 2>/dev/null) || echo FAIL $v; done
(cd w/tr && g++ $F -DDOUBLE -c main.cpp -o main.o 2>/dev/null && g++ $F -DDOUBLE main.o -o solve && g++ $F -DDOUBLE -S main.cpp -o main.s 2>/dev/null) || echo FAIL tr
for v in ours int i64 tr; do objdump -d --no-show-raw-insn -C w/$v/solve > w/$v/solve.dis; done
for inst in gr202 lin318; do
  cp $D/dm_$inst $D/w/distance_matrix
  for r in 1 2 3; do
    for v in ours int i64 tr; do
      out=$(cd w/$v && ./solve)
      echo "T $inst $v $r $(echo "$out" | grep -E 'wall clock' | awk '{print $NF}') $(echo "$out" | grep -o 'COST: *[0-9]*')"
    done
  done
done
echo FIM
