#!/bin/bash

#PBS -N grubert_hi-c-run
#PBS -o grubert_hi-c-output
#PBS -e grubert_hi-c-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=12
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

cd /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/godmc_phase2_analysis/13_hi-c
Rscript grubert_hic_run.r

####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.