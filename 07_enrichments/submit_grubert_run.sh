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

Rscript grubert_hic_new.r

####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.