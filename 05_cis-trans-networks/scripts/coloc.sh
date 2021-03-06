#!/bin/bash

#SBATCH --job-name=coloc
#SBATCH --nodes=1 --mem=2G --time=0-04:00:00
# #SBATCH --array=1-930
#SBATCH --array=28,629
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --partition=veryshort

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript coloc.r ${i} 1860

