#!/bin/bash

#PBS -N reformatAnnotation
#PBS -o /panfs/panasas01/shared-godmc/job_report/reformatannotation.o
#PBS -e /panfs/panasas01/shared-godmc/job_report/reformatannotation.e
#PBS -l walltime=100:00:00
# PBS -t 1
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash
# PBS -q himem

echo "Running on ${HOSTNAME}"

cd ~/repo/godmc_phase2_analysis/09_garfield/scripts
Rscript reformatAsAnnotation_Bristol.r

