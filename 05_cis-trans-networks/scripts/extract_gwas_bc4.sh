#!/bin/bash

#SBATCH --job-name=ml
#SBATCH --array=0-135
#SBATCH --nodes=1 --cpus-per-task=1 --time=0-00:30:00
#SBATCH --partition=mrcieu
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --mem=40G

echo "Running on ${HOSTNAME}"
module add languages/r/3.4.4

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}
# i=$((i + 1000))

# cd ${HOME}/mr-eve/gwas-instrument-subsets/scripts
ids=($(cat ../data/mrbids.txt))

echo ${#ids[@]}
id=`echo ${ids[$i]}`
echo $id
if [ -z "$id" ]
then
	echo "outside range"
	exit
fi

./extract_gwas.sh $id


