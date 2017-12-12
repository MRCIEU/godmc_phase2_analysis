#!/bin/bash

#SBATCH --job-name=extr
#SBATCH --nodes=1 --mem=20G --time=0-04:00:00
#SBATCH --partition=veryshort

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

cd ~/godmc/godmc_phase2_analysis/05_cis-trans-networks/scripts

Rscript extract_gwas.r

