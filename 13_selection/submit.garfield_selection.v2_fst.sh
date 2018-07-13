#!/bin/bash
  
#PBS -N garfield-fst   
#PBS -o garfield-fst-output
#PBS -e garfield-fst-error
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
echo fst
./garfield_fst

echo fst_cis
./garfield_fst_cis

echo fst_trans
./garfield_fst_trans

echo fst_ambivalent
./garfield_fst_ambivalent

