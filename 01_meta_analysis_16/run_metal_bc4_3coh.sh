#!/bin/bash

#SBATCH --job-name=meta_16_3cohorts
#SBATCH --array=1-962
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-00:30:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=100G

set -e

# Need to gunzip the results_16.tgz for each cohort
# change ARIES name to 00ARIES so that it processes first
# Make sure Dunedin38 excluded - it is a duplicate of Dunedin26
# Takes about 2 minutes to create each GWAMA file for each cohort


echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}
cd /mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/01_meta_analysis_16

metal="metal_bc4"

cohort_dir="/mnt/storage/home/epzjlm/repo/godmc_phase2_analysis/data/16_cohorts/"
metal_dir="../scratch_3coh/16_${i}"
result_dir="../results_3coh/16_3cohorts"

metal_in="16_${i}.in"

mkdir -p ${result_dir}
mkdir -p ${metal_dir}/scratch
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
	cohortname="scratch/$(basename "$cohort" .tar)"
	echo $cohortname
	mkdir -p ${cohortname}_${i}
	cd ${cohortname}_${i}
	tar xvf ${cohort} results/16/results_${i}.gz
	cd ../../
	mv ${cohortname}_${i}/results/16/results_${i}.gz ${metal_dir}/${cohortname}_${i}.gz
	rm -r ${cohortname}_${i}

	echo "PROCESS ${cohortname}_${i}.gz" >> ${metal_dir}/${metal_in}

done

echo "" >> ${metal_dir}/${metal_in}
echo "ANALYZE RANDOM" >> ${metal_dir}/${metal_in}
echo "CLEAR" >> ${metal_dir}/${metal_in}
echo "QUIT" >> ${metal_dir}/${metal_in}

cp ${metal} ${metal_dir}
cd ${metal_dir}

./${metal} ${metal_in} > 16_${i}.txt.log
mv 16_${i}1.txt 16_${i}.txt
mv 16_${i}1.txt.info 16_${i}.txt.info
gzip 16_${i}.txt
gzip 16_${i}.txt.log


## Remove bad alleles

zgrep "WARNING: Bad alleles for marker" 16_${i}.txt.log.gz | cut -d " " -f 7 | sed -E "s/'|,//g" | cut -d "_" -f 1 | sort -u > 16_${i}.txt.badlist
nbad=`cat 16_${i}.txt.badlist | wc -l`
echo "Removing ${nbad} bad SNPs"
zfgrep -vf 16_${i}.txt.badlist 16_${i}.txt.gz | gzip -c > 16_${i}.txt.gz2
mv 16_${i}.txt.gz2 16_${i}.txt.gz

cd -
mv ${metal_dir}/16_${i}.txt.* ${result_dir}
rm -r ${metal_dir}

# GWAMA --filelist ${metal_file} --quantitative

