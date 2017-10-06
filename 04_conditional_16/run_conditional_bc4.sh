#!/bin/bash

#SBATCH --job-name=cond_16
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --time=0-24:00:00
#SBATCH --array=1-962%100
#SBATCH --output=job_reports/slurm-%A_%a.out



set -e

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript conditional.r ${i} ../results/16/16_${i}_conditional.rdata

####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


