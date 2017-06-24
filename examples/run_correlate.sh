#!/bin/bash
# This script generates profile given the receptor and 
# the ligand. Then it scores the conformations provided
# in ft the files. Takes 5-20 min.

map=../prms/pdb_formfactor_mapping_clean.prm
prm=../prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec

lmax=15
c1=1.0
c2=1.0

pdbdir=others_example

rec=${pdbdir}/1a2k_r_u_nmin.pdb
lig=${pdbdir}/1a2k_l_u_nmin.pdb

rm_file=${pdbdir}/rot70k.0.0.6.jm.prm
ft_file1=${pdbdir}/ft.000.00
ft_file2=${pdbdir}/ft.001.00
ft_file3=${pdbdir}/ft.002.00

outdir=correlate_result
mkdir -p ${outdir}

chi=${outdir}/chi_scores
exp=${outdir}/ref_saxs_profile

# Running

len1=$(wc -l $ft_file1 | awk '{print $1}')
len2=$(wc -l $ft_file2 | awk '{print $1}')
len3=$(wc -l $ft_file3 | awk '{print $1}')

# Create reference SAXS profile
echo -ne "\n\e[93mGenerate reference profile\033[0m\n"
../build/tools/single_saxs $map $prm $rec $lig $c1 $c2 $lmax $exp
echo -ne "\e[93mProfile generated \033[0m\n"

# Score all conformations in ft files provided the reference profile
ft_file=${outdir}/ft_combo
cat $ft_file1 $ft_file2 $ft_file3 > $ft_file

echo -ne "\n\e[93mScore conformations with FFT-SAXS\033[0m\n"
mpirun -np 4 ../build/tools/correlate $map $prm $ft_file $rm_file $rec $lig $exp $lmax $chi

rm -f ${ft_file}

sort -k1 -n ${chi} -o ${chi}
tail -n${len3} ${chi} | cut -f2,3,4,5 > ${outdir}/chi.002.00
head -n${len1} ${chi} | cut -f2,3,4,5 > ${outdir}/chi.000.00
tail -n$((len2 + len3)) ${chi} | head -n${len2} | cut -f2,3,4,5 > ${outdir}/chi.001.00

echo -ne "\n\e[93mDone \033[0m\n"
