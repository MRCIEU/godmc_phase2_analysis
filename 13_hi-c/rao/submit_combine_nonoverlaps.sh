#!/bin/bash

#PBS -N nonoverlaps_combine
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/nonoverlaps_combine-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/nonoverlaps_combine-error
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

for i in {1..1000}; do
    cd /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/data/nonoverlaps/perm$i
    rm -f oe_data* bait_data* all_data_perm* nodups_data_perm*
    Rscript /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/combine_nonoverlaps.r $i
done
####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
