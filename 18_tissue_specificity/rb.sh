#!/bin/bash

#PBS -N rb
#PBS -o rb.o
#PBS -e rb.e

# PBS -t 1-23
#PBS -S /bin/bash
#! Resources requested:
#PBS -l nodes=1:ppn=8,walltime=12:00:00

#mydir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/01_meta_analysis_16"
mydir="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/18_tissue_specificity"
cd $mydir

#R CMD BATCH tissue.R tissue.Rout
R CMD BATCH calculate_rb.R calculate_rb.Rout

