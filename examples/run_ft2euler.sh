#!/bin/bash
# This script generates list of dimer conformations in Euler coordinates
# out of ft file.

pdbdir=dimer1

rec=${pdbdir}/r_native.pdb
lig=${pdbdir}/l_moved.pdb

rm_file=${pdbdir}/rot70k.0.0.6.jm.prm
ft_file=${pdbdir}/ft.000.00

euler=euler_list

echo -ne "\e[93mGenerate Euler list\033[0m\n"
../build/tools/ft2euler $ft_file $rm_file $rec $lig $euler
echo -ne "\n\e[93mDone \033[0m\n"
