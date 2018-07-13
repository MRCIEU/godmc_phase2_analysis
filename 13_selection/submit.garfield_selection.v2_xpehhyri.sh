#!/bin/bash
  
#PBS -N garfield-xpehhyri   
#PBS -o garfield-output-xpehhyri
#PBS -e garfield-error-xpehhyri
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
#echo xpehhyri
#./garfield_xpehhyri
echo xpehhyri_trans
./garfield_xpehhyri_trans
echo xpehhyri_cis
./garfield_xpehhyri_cis
echo xpehhyri_ambivalent
./garfield_xpehhyri_ambivalent

