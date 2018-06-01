#!/bin/bash
  
#PBS -N garfield-ihs   
#PBS -o garfield-ihs-output
#PBS -e garfield-ihs-error
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
#echo ihs
#./garfield_ihs
echo ihs_cis
./garfield_ihs_cis
echo ihs_trans
./garfield_ihs_trans
echo ihs_ambivalent
./garfield_ihs_ambivalent


