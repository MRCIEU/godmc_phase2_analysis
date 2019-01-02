#!/bin/bash

#SBATCH --job-name=enr
#SBATCH --time=0-10:30:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A.out
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28

echo "Running on ${HOSTNAME}"
module add languages/r/3.4.4

Rscript gwas_enrichment.r


