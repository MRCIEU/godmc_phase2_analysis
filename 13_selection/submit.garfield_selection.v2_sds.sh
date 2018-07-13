#!/bin/bash
  
#PBS -N garfield-sds   
#PBS -o garfield-sds-output
#PBS -e garfield-sds-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
# PBS -t 1-22

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2
echo sds
./garfield_sds
echo sds_cis
./garfield_sds_cis
echo sds_trans
./garfield_sds_trans
echo sds_ambivalent
./garfield_sds_ambivalent

