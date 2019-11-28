#!/bin/bash

#PBS -N mstat
#PBS -o job_report/mstat.o
#PBS -e job_report/mstat.e

# PBS -t 1-23
#PBS -S /bin/bash
#! Resources requested:
#PBS -l nodes=1:ppn=8,walltime=10:00:00

#mydir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/01_meta_analysis_16"
mydir="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/18_tissue_specificity"
cd $mydir

#R CMD BATCH mstat.control.R mstat.Rout

#R CMD BATCH extractclumpedsnps.R extractclumpedsnps.Rout 
#plink --bfile /panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered --extract clumpedsnps.txt --indep 50 5 1.010101 --out indep.clump
#R CMD BATCH extractchunks.R extractchunks.Rout
#sh extractchr20chunks.sh

R CMD BATCH heterogeneity.R heterogeneity.Rout
mv *tif ./images
