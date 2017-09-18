#!/bin/bash

#SBATCH --job-name=lola
#SBATCH --nodes=1 --mem=60G --ntasks=10 --time=0-20:00:00
#SBATCH --output=job_reports/slurm-%A.out

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

Rscript lola.r

