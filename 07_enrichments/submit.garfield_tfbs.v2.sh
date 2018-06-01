#!/bin/bash
  
#PBS -N garfield-tfbs   
#PBS -o garfield-tfbs-output
#PBS -e garfield-tfbs-error
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=6
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
perl -pe 's/T lymphocyte/T_lymphocyte/g' <../garfield-data/annotation-tfbs/link_file_tfbs.txt |perl -pe 's/B lymphocyte/B_lymphocyte/g'|perl -pe 's/ cell/\_cell/g' |perl -pe 's/ carcinoma/\_carcinoma/g'|perl -pe 's/Blood vessel/Blood_vessel/g' |perl -pe 's/Spinal cord/Spinal_cord/g'| perl -pe 's/Choroid plexus/Choroid_plexus/g' >../garfield-data/annotation-tfbs/link_file_tfbs2.txt

./garfield_tfbs



