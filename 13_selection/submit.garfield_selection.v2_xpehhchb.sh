#!/bin/bash
  
#PBS -N garfield-chb   
#PBS -o garfield-chb-output
#PBS -e garfield-chb-error
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
#echo xpehhchb
#./garfield_xpehhchb
echo xpehhchb_trans
./garfield_xpehhchb_trans
echo xpehhchb_cis
./garfield_xpehhchb_cis
echo xpehhchb_ambivalent
./garfield_xpehhchb_ambivalent
