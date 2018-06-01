#!/bin/bash
  
#PBS -N cpgmean-tfbs  
#PBS -o cpgmeantfbs-output
#PBS -e cpgmeantfbs-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash
# PBS -t 1-23
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
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments
R CMD BATCH --no-save --no-restore cpgsbytfbs.R cpgsbytfbs.Rout



