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

declare -a cohorts=(
"00ARIES"
"ALS_Batch1"
"ALS_Batch2"
"BAMSE"
"BASICMAR"
"BSGS"
"BorninBradford"
"DunedinAge26"
"DunedinAge38"
"EPICOR"
"Erisk"
"FTC"
"GLAKU"
"GOYA"
"GSK"
"GenerationR_Study"
"INMA"
"IOW3g"
"InterAct"
"LBC21"
"LBC36"
"Leiden_Longevity_Study"
"MARTHA"
"NTR"
"PIAMA"
"PRECISESADS"
"PREDO"
"Phase1_SCZ"
"Phase2_SCZ"
"Project_MinE_s27"
"Raine"
"Rotterdam_Study"
"SABRE"
"TwinsUK"
"UnderstandingSociety"
"as_cc"
"ccg"
)


echo "Running on ${HOSTNAME}"

i=${PBS_ARRAYID}

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

cohort=${cohorts[${i}]}
echo ${cohort}

cd ../
rtdr=`pwd`
cd -
new_cohort_dir="${rtdr}/data/16"
mkdir -p ${new_cohort_dir}
flipfiles="${rtdr}/data/flipped_snps"

# for cohort in "${cohorts[@]}"
# do
echo $cohort
scratch_dir="${rtdr}/scratch/16_${cohort}"
rm -rf ${scratch_dir}
mkdir -p ${scratch_dir}
cd ${scratch_dir}
tar xf ${rtdr}/data/16_raw/${cohort}_16.tar
# for i in {1..962}
# do
# 	echo ${i}
# 	Rscript ${rtdr}/01_meta_analysis_16/flip_alleles.r results/16/results_${i}.gz ${flipfiles}/${cohort}_data.easyqc.flipped.SNPs.txt ${rtdr}/data/ref/eur2.bim.rdata
# done
parallel -j12 Rscript ${rtdr}/01_meta_analysis_16/flip_alleles.r results/16/results_{}.gz ${flipfiles}/${cohort}_data.easyqc.flipped.SNPs.txt ${rtdr}/data/ref/eur2.bim.rdata ::: {1..962}

ls -lrt results/16/results* | head
tar cf ${new_cohort_dir}/${cohort}_16.tar *
cd ${rtdr}/01_meta_analysis_16
rm -rf ${scratch_dir}



