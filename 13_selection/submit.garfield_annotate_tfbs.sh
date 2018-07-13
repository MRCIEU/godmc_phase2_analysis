#!/bin/bash
  
#PBS -N garfield-annotation   
#PBS -o garfield-annotation-output
#PBS -e garfield-annotation-error
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash
# PBS -t 1-25

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
./garfield_annotate_uk10k.sh
