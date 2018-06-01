#!/bin/bash
  
#PBS -N garfield-segmentations-cis   
#PBS -o garfield-segmentations-cis-output
#PBS -e garfield-segmentations-cis-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
#PBS -t 1-25

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
perl -pe 's/STATE=1/STATE='$i'/g' < garfield_segmentations_cis > garfield_segmentations_cis${i}
chmod +x garfield_segmentations_cis${i}
./garfield_segmentations_cis${i}
