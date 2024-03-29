#!/bin/bash
  
#PBS -N garfield-segmentations-cis-sds   
#PBS -o garfield-segmentations-cis-sds-output
#PBS -e garfield-segmentations-cis-sds-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=2
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
perl -pe 's/STATE=1/STATE='$i'/g' < garfield_segmentations_cis_sds > garfield_segmentations_cis_sds${i}
chmod +x garfield_segmentations_cis_sds${i}
./garfield_segmentations_cis_sds${i}
