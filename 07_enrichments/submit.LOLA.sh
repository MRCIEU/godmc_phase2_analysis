#!/bin/bash
  
#PBS -N LOLA   
#PBS -o LOLA-output
#PBS -e LOLA-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
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

cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/

R CMD BATCH --no-save --no-restore cpg_enrichments_lola.r cpg_enrichments_lola.Rout

R CMD BATCH --no-save --no-restore cpg_enrichments_lola_ext.r cpg_enrichments_lola_ext.Rout

R CMD BATCH --no-save --no-restore cpg_enrichments_lola_chromstates.r cpg_enrichments_lola_chromstates.Rout


