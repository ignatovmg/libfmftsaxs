#!/bin/bash
# This script generates SAXS profile
# with specified parameters c1 and c2

map=../prms/pdb_formfactor_mapping_clean.prm
prm=../prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec

rec=dimer1/r_native.pdb
lig=dimer1/l_native.pdb

c1=1.0 # range [0.96, 1.04]
c2=0.0 # range [-2, 4]
lmax=15

profile=dimer1_saxs_profile

echo -ne "\e[93mGenerate profile with c1 = ${c1}, c2 = ${c2} \033[0m\n"
../build/tools/single_saxs $map $prm $rec $lig $c1 $c2 $lmax $profile

echo -ne "\n\e[93mDone \033[0m\n"
