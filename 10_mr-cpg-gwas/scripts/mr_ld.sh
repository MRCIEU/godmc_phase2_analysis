#!/bin/bash

#SBATCH --job-name=mr_ld
#SBATCH --nodes=1 --mem=10G --time=0-40:00:00
#SBATCH --array=301-701
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=mrcieu

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript mr_ld.r ${i}

