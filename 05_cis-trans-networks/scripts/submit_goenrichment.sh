#!/bin/bash
  
#PBS -N goenrich  
#PBS -o goenrich-output
#PBS -e goenrich-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash
# PBS -t 1-23
# PBS -q himem

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/05_cis-trans-networks/scripts
R CMD BATCH --no-save --no-restore KEGG_gsa_enrichment.R KEGG_gsa_enrichment.Rout

