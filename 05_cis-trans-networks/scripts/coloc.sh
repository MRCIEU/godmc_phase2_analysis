#!/bin/bash

#SBATCH --job-name=coloc
#SBATCH --nodes=1 --mem=4G --time=0-20:00:00
#SBATCH --array=203,304,65,66,67,68,69,70,71,72,73,74,75,76

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript coloc.r ${i} 100

