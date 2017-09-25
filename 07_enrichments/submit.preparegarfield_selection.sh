#!/bin/bash
  
#PBS -N prepare-garfield   
#PBS -o garfieldprep-output
#PBS -e garfieldprep-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
#PBS -t 1-22

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd ~/repo/godmc_phase2_analysis/07_enrichments

R CMD BATCH --no-save --no-restore '--args '$i'' preparegarfield_selection.R preparegarfield_selection$i.Rout

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield
./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl_ihs/chr$i ../garfield-data/annotation_selection/chr$i -1 > ../garfield-data/prep_mqtl_ihs/chr$i
./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/transmqtl_ihs/chr$i ../garfield-data/annotation_selection_trans/chr$i -1 > ../garfield-data/prep_transmqtl_ihs/chr$i

./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl_fst/chr$i ../garfield-data/annotation_selection/chr$i -1 > ../garfield-data/prep_mqtl_fst/chr$i
./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/transmqtl_fst/chr$i ../garfield-data/annotation_selection_trans/chr$i -1 > ../garfield-data/prep_transmqtl_fst/chr$i

./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl_xpehhchb/chr$i ../garfield-data/annotation_selection/chr$i -1 > ../garfield-data/prep_mqtl_xpehhchb/chr$i
./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/transmqtl_xpehhchb/chr$i ../garfield-data/annotation_selection_trans/chr$i -1 > ../garfield-data/prep_transmqtl_xpehhchb/chr$i

./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl_xpehhyri/chr$i ../garfield-data/annotation_selection/chr$i -1 > ../garfield-data/prep_mqtl_xpehhyri/chr$i
./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/transmqtl_xpehhyri/chr$i ../garfield-data/annotation_selection_trans/chr$i -1 > ../garfield-data/prep_transmqtl_xpehhyri/chr$i

./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/mqtl_sds/chr$i ../garfield-data/annotation_selection/chr$i -1 > ../garfield-data/prep_mqtl_sds/chr$i
./garfield-prep ../garfield-data/tags/r01/chr$i ../garfield-data/tags/r08/chr$i ../garfield-data/maftssd/chr$i ../garfield-data/pval/transmqtl_sds/chr$i ../garfield-data/annotation_selection_trans/chr$i -1 > ../garfield-data/prep_transmqtl_sds/chr$i


#BEFORE RUNNING THIS SCRIPT:
# you need to run 
# preparegarfield_selection0.R
# compile.snpcpgpval.sh
