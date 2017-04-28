#!/bin/bash

#PBS -N meta_16
#PBS -o job_reports/meta_16-output
#PBS -e job_reports/meta_16-error
#PBS -t 901-962
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash



set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd ~/repo/godmc_phase2_analysis/03_clumping_16
Rscript clump.r ${i} ../results/16/16_${i}_clumped.rdata

