#!/bin/bash

#SBATCH --job-name=flip_16
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --time=0-04:00:00
#SBATCH --array=10,20
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --queue=veryshort

set -e


echo "Running on ${HOSTNAME}"

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

cohort="EPICOR"
echo ${cohort}

cd ../
rtdr=`pwd`
cd -
renamed_orig="${rtdr}/data/16_raw/${cohort}_16.tar.orig"
scratch_dir="${rtdr}/scratch/16_${cohort}"
rm -rf ${scratch_dir}
mkdir -p ${scratch_dir}

# rename SNPs in indel list
wc -l ../data/16_raw/EPICOR_INDEL_list.txt
# sed -E 's/(:I)|(:D)/:INDEL/g' ../data/16_raw/EPICOR_INDEL_list.txt > ../data/16_raw/EPICOR_INDEL_list.txt2

cd ${scratch_dir}
tar xf ${rtdr}/data/16_raw/${cohort}_16.tar
# for i in {1..962}
# do
# 	echo ${i}
# 	Rscript ${rtdr}/01_meta_analysis_16/fix_epicor.r results/16/results_${i}.gz ${rtdr}/data/16_raw/EPICOR_INDEL_list.txt
# done
parallel -j12 Rscript ${rtdr}/01_meta_analysis_16/fix_epicor.r results/16/results_{}.gz ${rtdr}/data/16_raw/EPICOR_INDEL_list.txt ::: {1..962}

ls -lrt results/16/results* | head

mv ${rtdr}/data/16_raw/${cohort}_16.tar ${renamed_orig}
tar cf ${rtdr}/data/16_raw/${cohort}_16.tar *
cd ${rtdr}/01_meta_analysis_16
rm -rf ${scratch_dir}


