#!/bin/bash
  
#PBS -N snpset   
#PBS -o snpsetgc-output
#PBS -e snpsetgc-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash
# PBS -t 1-23
# PBS -t 1
#PBS -q himem

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
#R CMD BATCH addcg_cpgfrq.R addcg_cpgfrq.Rout
#R CMD BATCH ld_regions_makecontrols.R ld_regions_makecontrols.Rout

#R CMD BATCH --no-save --no-restore '--args '$i'' ld_regions_makecontrols.R ld_regions_makecontrols.Rout 
R CMD BATCH --no-save --no-restore '--args '$i'' addcg_cpgfrq.R addcg_cpgfrq$i.Rout 
#R CMD BATCH --no-save --no-restore preparegarfield_selection0.R preparegarfield_selection0.Rout
