#!/bin/bash
  
#PBS -N garfield-sds   
#PBS -o garfield-sds-output
#PBS -e garfield-sds-error
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
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

#R CMD BATCH --no-save --no-restore '--args '$i'' preparegarfield_epigeneticannotation.R preparegarfield$i.Rout

cd /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2
./garfield_sds_epigenetic


