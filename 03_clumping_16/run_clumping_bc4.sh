#!/bin/bash

#SBATCH --job-name=clump_16
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --time=0-12:00:00
#SBATCH --array=1-962%100
#SBATCH --output=job_reports/slurm-%A_%a.out

set -e

echo "Running on ${HOSTNAME}"

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

# cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/03_clumping_16
Rscript clump_bc4.r \
	${i} \
	../results/16/16_${i}_clumped.rdata \
	1e-4 \
	5e-8 \
	0.0001 \
	5000 \
	1000000
