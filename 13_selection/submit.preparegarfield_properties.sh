#!/bin/bash
  
#PBS -N prepare-garfield   
#PBS -o garfieldprep-output
#PBS -e garfieldprep-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
#PBS -t 23 

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

R CMD BATCH --no-save --no-restore '--args '$i'' cg_cpgfrq_GARFIELD.R cg_cpgfrq_GARFIELD$i.Rout



#BEFORE RUNNING THIS SCRIPT:
# you need to run 
# preparegarfield_selection0.R
# compile.snpcpgpval.sh

dir="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/maftssd/"
cd $dir

for i in `seq 1 23`; do
echo $i
##mv chr$i chr$i.orig
awk '{print $1, $2,$3,$4,$5,$6}' <chr$i.tmp >chr$i
done