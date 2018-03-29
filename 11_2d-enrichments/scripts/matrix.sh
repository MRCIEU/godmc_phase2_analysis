#!/bin/bash

#SBATCH --job-name=2d
#SBATCH --nodes=1 
#SBATCH --mem=60G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=0-20:00:00
#SBATCH --array=0-100%50
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript matrix.r ${i}

