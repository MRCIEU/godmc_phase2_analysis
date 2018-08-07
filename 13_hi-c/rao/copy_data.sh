#!/bin/bash

#PBS -N copy_data2
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/copy_data2-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/copy_data2-error
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
## PBS -t 253
#PBS -t 201-253
## PBS -t 201-253
#PBS -S /bin/bash

# To Do. 
#- copy data files to data directory

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

#1. copy data NORM files
#echo ${x}.NORM
#wc -l ${x}.NORM
#cp ${x}.NORM ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/norm_files/
#wc -l ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/norm_files/${x}.NORM

#2. copy data overlap files
#echo data_${x}.Rdata
#cp data_${x}.Rdata ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/overlaps/

#3. copy permutation  nonoverlap files
for i in {1..1000}; do
echo $i
echo nondata_${x}_perm_$i.Rdata
cp -f nondata_${x}_perm_$i.Rdata ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/nonoverlaps/perm$i

echo "$x now synced"
done
done