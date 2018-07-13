#!/bin/bash
  
#PBS -N prepare-tfbs  
#PBS -o garfieldtfbs-output
#PBS -e garfieldtfbs-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash
#PBS -t 1-23

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments
R CMD BATCH --no-save --no-restore '--args '$i'' formatbedforGARFIELD.R formatbedforGARFIELD$i.Rout



