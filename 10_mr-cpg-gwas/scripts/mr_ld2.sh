#!/bin/bash

#SBATCH --job-name=mr_ld2
#SBATCH --nodes=1 --mem=2G --ntasks=1 --time=0-12:00:00
#SBATCH --array=1-999
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript mr_ld2.r ${i}

