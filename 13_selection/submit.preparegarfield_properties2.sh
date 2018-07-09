#!/bin/bash
  
#PBS -N prepare-garfield   
#PBS -o garfieldprep-output
#PBS -e garfieldprep-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
# PBS -t 1-23

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

#R CMD BATCH --no-save --no-restore addcg_cpgfrq.R addcg_cpgfrq.Rout
R CMD BATCH --no-save --no-restore preparegarfield_selection0.R 


#BEFORE RUNNING THIS SCRIPT:
# you need to run 
# preparegarfield_selection0.R
# compile.snpcpgpval.sh
