#!/bin/bash

#PBS -N meta_16
#PBS -o job_reports/meta_16-output
#PBS -e job_reports/meta_16-error
#PBS -t 1-100
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash



set -e

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/04_conditional_16

Rscript conditional.r ${i} ../results/16/16_${i}_conditional.rdata



####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


