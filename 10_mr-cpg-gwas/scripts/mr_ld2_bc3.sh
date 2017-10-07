#!/bin/bash

#PBS -N mr_ld2
#PBS -o job_reports/mr_ld2-output
#PBS -e job_reports/mr_ld2-error
#PBS -t 1-7010
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash



set -e

# Need to gunzip the results_16.tgz for each cohort
# Takes about 2 minutes to create each GWAMA file for each cohort


echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/10_mr-cpg-gwas/scripts/

Rscript mr_ld2.r ${i}

