#!/bin/bash

#SBATCH --job-name=coloc
#SBATCH --nodes=1 --mem=20G --time=0-20:00:00
#SBATCH --array=1-158
# #SBATCH --array=7,29,88,99,103,154,158
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript coloc.r ${i} 100

