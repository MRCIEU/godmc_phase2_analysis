#!/bin/bash

#PBS -N mstat
#PBS -o job_report/mstat.o
#PBS -e job_report/mstat.e

# PBS -t 1-23
#PBS -S /bin/bash
#! Resources requested:
#PBS -l nodes=1:ppn=8,walltime=10:00:00

#mydir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/01_meta_analysis_16"
mydir="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/01_meta_analysis_16"
cd $mydir

#R CMD BATCH mstat.control.R mstat.Rout

#sh extractchr20chunks.sh
R CMD BATCH mstat.R mstat.Rout
