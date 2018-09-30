#!/bin/bash

#PBS -N overlaps
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/overlaps-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/overlaps-error
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -t 1-100
#PBS -S /bin/bash

#1-253 chr_list

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

IFS=_ read var1 var2 <<< "$x"
echo $var1
echo $var2

#rm -f hi_c_ranges_${x}.rdata
#rm -f trans_mqtl_granges_${x}.rdata

Rscript /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/rao_hic_overlaps.r $x $var1 $var2
done
####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

