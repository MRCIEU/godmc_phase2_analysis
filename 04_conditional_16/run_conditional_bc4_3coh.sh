#!/bin/bash
#SBATCH --job-name=cond_16_3cohorts
#SBATCH --array=1-962
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-4:00:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=4G

set -e

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}
cd /mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/04_conditional_16
Rscript /mnt/storage/home/epzjlm/repo/godmc_phase2_analysis/04_conditional_16/conditional_3coh.r ${i} ../results_3coh/16_3cohorts/16_${i}_conditional.rdata

####

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.


