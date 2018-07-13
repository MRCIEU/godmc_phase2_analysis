#!/bin/bash
  
#PBS -N snpset1   
#PBS -o snpset1-output
#PBS -e snpset1-error
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
# PBS -t 1-10

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd ~/repo/godmc_phase2_analysis/13_selection

R CMD BATCH --no-save --no-restore '--args '$i'' ld_regions_makecontrols.R ld_regions_makecontrols.Rout 
