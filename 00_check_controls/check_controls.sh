#!/bin/bash

#PBS -N controls
#PBS -o job_reports/controls-output
#PBS -e job_reports/controls-error
#PBS -t 1-962
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash



set -e

# Need to gunzip the results_16.tgz for each cohort
# Takes about 2 minutes to create each GWAMA file for each cohort


echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}



cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/00_check_controls

cohort_dir="../data/16/"
result_dir="../results/controls"

mkdir -p ${result_dir}

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	mkdir -p ${cohortname}
	cd ${cohortname}
	tar xvf ../${cohort} results/16/logs_b/log.txt
	tar xvf ../${cohort} results/16/control/cg07959070_qqplot.png
	tar xvf ../${cohort} results/16/control/cg07959070_manhattan.png
	tar xvf ../${cohort} results/16/control/cg07959070_without_chr22_qqplot.png

	mv results/16/logs_b/log.txt ../${result_dir}/${cohortname}_log.txt
	mv results/16/control/cg07959070_qqplot.png ../${result_dir}/${cohortname}_cg07959070_qqplot.png
	mv results/16/control/cg07959070_manhattan.png ../${result_dir}/${cohortname}_cg07959070_manhattan.png
	mv results/16/control/cg07959070_without_chr22_qqplot.png ../${result_dir}/${cohortname}_cg07959070_without_chr22_qqplot.png
	cd ../
	rm -r ${cohortname}

	echo ${cohortname}
	grep "lambda value for GWAS" ${result_dir}/${cohortname}_log.txt | cut -d " " -f 5
done

