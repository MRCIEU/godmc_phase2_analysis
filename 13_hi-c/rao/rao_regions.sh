#!/bin/bash

#PBS -N regions-Rao
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/regions-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/regions-error
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -t 1-100
#PBS -S /bin/bash

# To Do. 
#- take contact point 1 + 1kb
#- take contact point 2 + 1kb

set -e

echo $PBS_ARRAYID

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

for NUM in $PBS_ARRAYID; do
dir=($(sed "${NUM}q;d" ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/chr_list.txt))
done

for x in ${dir[@]}; do
echo $x
cd ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/$x/MAPQGE30/
pwd

IFS=_ read var1 var2 <<< "$x"
echo $var1
echo $var2

#1. add 1kb regions
rm -f ${x}.NORM
awk -v a="$var1" -v b="$var2" '{ print $1, "\t", $1+1000, "\t", $2, "\t", $2+1000, "\t", $3, "\t", a, "\t", b }' ${x}_1kb.NORMobserved > ${x}.NORM
done



