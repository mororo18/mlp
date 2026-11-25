#!/usr/bin/env bash
# Perfil do LuaJIT (-jp) no gr202: nossa com LJ1 e tradução com L1 e com L1-L5.
cd "$(dirname "$0")"; D=$PWD
cp $D/dm_gr202 $D/w/distance_matrix
for d in cpp_lua.luajit_subseq_array tr_all_hoist_carry; do
  (cd $D/w/$d && luajit -jp=Fl2i1 main.lua > $D/jp_gr202_$d.txt 2>&1)
  (cd $D/w/$d && luajit -jp=v main.lua > $D/jpv_gr202_$d.txt 2>&1)
done
echo FIM
