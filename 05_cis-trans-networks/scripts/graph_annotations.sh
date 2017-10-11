#!/bin/bash

#SBATCH --job-name=annot
#SBATCH --nodes=1 --mem=4G --time=0-20:00:00
#SBATCH --array=1-57
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript graph_annotations.r ${i}

