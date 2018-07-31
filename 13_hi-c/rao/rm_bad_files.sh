#!/bin/bash
#PBS -N rmBadFiles
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/lscount2-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/lscount2-error
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -t 1-100
#PBS -S /bin/bash

set -e

echo $PBS_ARRAYID

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

for NUM in $PBS_ARRAYID; do
dir=($(sed "${NUM}q;d" ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/chr_list.txt))
done

for x in ${dir[@]}; do
echo $x
cd ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/$x/MAPQGE30
#rm -rf non*
#rm -rf *chr18_chrX*
ls -l *perm* | wc -l
#cp -rv /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/$x /projects/MRC-IEU/research/data/godmc/_devs/GODMC_Analysis/data/Hi-C/rao
#ls -l
done

