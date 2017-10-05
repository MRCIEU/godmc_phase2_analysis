#!/bin/bash

#SBATCH --job-name=meta_16
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00
#SBATCH --array=11-962%100
#SBATCH --output=job_reports/slurm-%A_%a.out


set -e

# Need to gunzip the results_16.tgz for each cohort
# change ARIES name to 00ARIES so that it processes first
# Make sure Dunedin38 excluded - it is a duplicate of Dunedin26
# Takes about 2 minutes to create each GWAMA file for each cohort


echo "Running on ${HOSTNAME}"

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

metal="metal_bc4"

cohort_dir="../data/16/"
metal_dir="../scratch/16_${i}"
result_dir="../results/16"
new_cohort_dir="../data/16_flipped"

metal_in="16_${i}.in"

mkdir -p ${result_dir}
mkdir -p ${new_cohort_dir}
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
	cohortname1="$(basename "$cohort" .tar)"
	mkdir -p ${new_cohort_dir}/${cohortname1}
	cohortname="scratch/$(basename "$cohort" .tar)"
	mkdir -p ${cohortname}_${i}
	cd ${cohortname}_${i}
	tar xvf ../../${cohort} results/16/results_${i}.gz
	cd ../../
	gunzip ${cohortname}_${i}/results/16/results_${i}.gz
	Rscript flip_alleles.r \
		${cohortname}_${i}/results/16/results_${i} \
		${new_cohort_dir}/${cohortname1}/results_${i}.gz
	rm ${cohortname}_${i}/results/16/results_${i}
	cp ${new_cohort_dir}/${cohortname1}/results_${i}.gz ${metal_dir}/${cohortname}_${i}.gz
	rm -r ${cohortname}_${i}

	echo "PROCESS ${cohortname}_${i}.gz" >> ${metal_dir}/${metal_in}

done

echo "" >> ${metal_dir}/${metal_in}
echo "ANALYZE RANDOM" >> ${metal_dir}/${metal_in}
echo "CLEAR" >> ${metal_dir}/${metal_in}
echo "QUIT" >> ${metal_dir}/${metal_in}

cp ${metal} ${metal_dir}
cd ${metal_dir}

./${metal} ${metal_in}
mv 16_${i}1.txt 16_${i}.txt
gzip 16_${i}.txt
cd -
mv ${metal_dir}/16_${i}.txt.* ${result_dir}
# rm -r ${metal_dir}

# GWAMA --filelist ${metal_file} --quantitative

