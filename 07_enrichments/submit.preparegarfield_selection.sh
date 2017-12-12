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

#datadir="/panfs/panasas01/shared-godmc/GARFIELD/garfield-data"
#cd /panfs/panasas01/shared-godmc/GARFIELDv2/garfield
#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/mqtl_ihs/chr$i $datadir/annotation_selection/chr$i -1 > $datadir/prep_mqtl_ihs/chr$i
#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/transmqtl_ihs/chr$i $datadir/annotation_selection_trans/chr$i -1 > $datadir/prep_transmqtl_ihs/chr$i

#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/mqtl_fst/chr$i $datadir/annotation_selection/chr$i -1 > $datadir/prep_mqtl_fst/chr$i
#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/transmqtl_fst/chr$i $datadir/annotation_selection_trans/chr$i -1 > $datadir/prep_transmqtl_fst/chr$i

#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/mqtl_xpehhchb/chr$i $datadir/annotation_selection/chr$i -1 > $datadir/prep_mqtl_xpehhchb/chr$i
#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/transmqtl_xpehhchb/chr$i $datadir/annotation_selection_trans/chr$i -1 > $datadir/prep_transmqtl_xpehhchb/chr$i

#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/mqtl_xpehhyri/chr$i $datadir/annotation_selection/chr$i -1 > $datadir/prep_mqtl_xpehhyri/chr$i
#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/transmqtl_xpehhyri/chr$i $datadir/annotation_selection_trans/chr$i -1 > $datadir/prep_transmqtl_xpehhyri/chr$i

#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/mqtl_sds/chr$i $datadir/annotation_selection/chr$i -1 > $datadir/prep_mqtl_sds/chr$i
#./garfield-prep $datadir/tags/r01/chr$i $datadir/tags/r08/chr$i $datadir/maftssd/chr$i $datadir/pval/transmqtl_sds/chr$i $datadir/annotation_selection_trans/chr$i -1 > $datadir/prep_transmqtl_sds/chr$i


#BEFORE RUNNING THIS SCRIPT:
# you need to run 
# preparegarfield_selection0.R
# compile.snpcpgpval.sh
