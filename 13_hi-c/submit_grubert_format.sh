#!/bin/bash

#PBS -N hi-c-clean
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Grubert2015/job_reports/clean-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Grubert2015/job_reports/clean-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=12
## PBS -t 0-99
#PBS -t 100-129
#PBS -S /bin/bash

set -e
printf -v nm %03d $PBS_ARRAYID

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

cd /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Grubert2015/

Rscript Grubert_Format.R $nm



####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.