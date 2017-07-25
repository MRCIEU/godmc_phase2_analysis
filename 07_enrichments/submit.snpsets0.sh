#!/bin/bash
  
#PBS -N snpset   
#PBS -o snpset-output
#PBS -e snpset-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=16
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

cd ~/repo/godmc_phase2_analysis/07_enrichments
R CMD BATCH ld_regions.r ld_regions.Rout
#R CMD BATCH ld_regions_makecontrols.R ld_regions_makecontrols.Rout

#R CMD BATCH --no-save --no-restore '--args '$i'' ld_regions_makecontrols.R ld_regions_makecontrols.Rout 
