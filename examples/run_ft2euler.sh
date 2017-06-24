#!/bin/bash
# This script generates list of dimer conformations in Euler coordinates
# out of ft file.

pdbdir=4g9s

rec=${pdbdir}/4g9s_r_native.pdb
lig=${pdbdir}/4g9s_l_moved.pdb

rm_file=${pdbdir}/rm.000.00
ft_file=${pdbdir}/ft.000.00

euler=euler_list

echo -ne "\e[93mGenerate Euler list\033[0m\n"
../build/tools/ft2euler $ft_file $rm_file $rec $lig $euler
echo -ne "\n\e[93mDone \033[0m\n"
