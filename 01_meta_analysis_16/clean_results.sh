#!/bin/bash

#SBATCH --job-name=meta_16
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --array=1-962
#SBATCH --time=0-02:00:00

set -e

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript clean_results.r ${i}
