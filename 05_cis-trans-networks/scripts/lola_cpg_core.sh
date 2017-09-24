#!/bin/bash

#SBATCH --job-name=lola_cpg_core
#SBATCH --nodes=1 --mem=120G --ntasks=5 --time=0-20:00:00
#SBATCH --output=job_reports/slurm-%A.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

Rscript lola_cpg_core.r

