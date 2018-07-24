#!/bin/bash

#SBATCH --job-name=difres
#SBATCH --nodes=1 
#SBATCH --mem=20G
#SBATCH --time=0-02:00:00
#SBATCH --array=0-1000
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --partition=mrcieu

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript difres.r ${i}

