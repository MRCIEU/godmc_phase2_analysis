#!/bin/bash

#SBATCH --job-name=gwas
#SBATCH --nodes=1 --mem=10G --time=0-40:00:00
#SBATCH --array=1-942
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript gwas.r ${i} 500 1e-5

