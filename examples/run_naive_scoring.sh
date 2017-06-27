#!/bin/bash
# This script scores conformations in ft file
# in straightforward one-by-one way. 
# Takes several hours.

map=../prms/pdb_formfactor_mapping_clean.prm
prm=../prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec

lmax=15
pdbdir=dimer1

rec=${pdbdir}/r_native.pdb
lig=${pdbdir}/l_moved.pdb

rm_file=${pdbdir}/rot70k.0.0.6.jm.prm
ft_file=${pdbdir}/ft.000.00

exp=${pdbdir}/ref_saxs # experimental profile

chi=chi_scores # output

echo -ne "\e[93mScore conformations straightforward\033[0m\n"
../build/tools/score_ft_naive $map $prm $ft_file $rm_file $rec $lig $exp $lmax $chi
echo -ne "\n\e[93mDone\033[0m\n"
