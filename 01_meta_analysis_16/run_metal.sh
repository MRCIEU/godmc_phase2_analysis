#!/bin/bash

#PBS -N meta_16
#PBS -o job_reports/meta_16-output
#PBS -e job_reports/meta_16-error
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



cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/01_meta_analysis_16

cohort_dir="../data/16/"
metal_dir="../scratch/16_${i}"
result_dir="../results/16"

metal_in="16_${i}.in"

mkdir -p ${result_dir}
mkdir -p ${metal_dir}
rm -f ${metal_dir}/${metal_in}
touch ${metal_dir}/${metal_in}



echo "OUTFILE 16_${i} .txt" >> ${metal_dir}/${metal_in}
echo "" >> ${metal_dir}/${metal_in}
echo "MARKER MARKERNAME" >> ${metal_dir}/${metal_in}
echo "ALLELE EA NEA" >> ${metal_dir}/${metal_in}
echo "EFFECT BETA" >> ${metal_dir}/${metal_in}
echo "FREQ EAF" >> ${metal_dir}/${metal_in}
echo "STDERR SE" >> ${metal_dir}/${metal_in}
echo "WEIGHT N" >> ${metal_dir}/${metal_in}
echo "CUSTOMVARIABLE TotalSampleSize" >> ${metal_dir}/${metal_in}
echo "LABEL TotalSampleSize as N" >> ${metal_dir}/${metal_in}
echo "SCHEME STDERR" >> ${metal_dir}/${metal_in}
echo "AVERAGEFREQ ON" >> ${metal_dir}/${metal_in}
echo "" >> ${metal_dir}/${metal_in}

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	mkdir -p ${cohortname}_${i}
	cd ${cohortname}_${i}
	tar xvf ../${cohort} results/16/results_${i}.gz
	cd ../
	mv ${cohortname}_${i}/results/16/results_${i}.gz ${metal_dir}/${cohortname}_${i}.gz
	rm -r ${cohortname}_${i}

	echo "PROCESS ${cohortname}_${i}.gz" >> ${metal_dir}/${metal_in}

done

echo "" >> ${metal_dir}/${metal_in}
echo "ANALYZE HETEROGENEITY" >> ${metal_dir}/${metal_in}
echo "CLEAR" >> ${metal_dir}/${metal_in}
echo "QUIT" >> ${metal_dir}/${metal_in}

cp metal ${metal_dir}
cd ${metal_dir}

./metal ${metal_in}
mv 16_${i}1.txt 16_${i}.txt
gzip 16_${i}.txt
cd -
mv ${metal_dir}/16_${i}.txt.* ${result_dir}
rm -r ${metal_dir}

# GWAMA --filelist ${metal_file} --quantitative

