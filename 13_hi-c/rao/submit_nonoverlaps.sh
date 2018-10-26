#!/bin/bash

#PBS -N nonoverlaps
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/nonoverlaps-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/nonoverlaps-error
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=2
#PBS -t 131-134,191,253,3,68,77,85
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

    for i in {1..1000}; do
	    Rscript --verbose /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/nonrao_hic_overlaps.r $x $var1 $var2 $i
    done
done
####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.