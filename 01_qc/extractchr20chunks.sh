#!bin/bash

results_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
cohort_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/16/"


cd $results_dir

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	echo $cohortname
	mkdir -p $cohort_dir/${cohortname}

		cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/01_meta_analysis_16
		file="chr20_chunks.txt"
		while read -r line
		do
			chunk="$line"
			echo "$chunk"

			cd $cohort_dir/${cohortname}
			tar -xvf ../${cohortname}.tar results/16/results_${chunk}.gz
		done < $file

###for i in `seq 1 25`; 
###do
###	tar -xvf ../${cohortname}.tar results/16/results_${i}.gz
###done
rm -r $results_dir/${cohortname}
mv $cohort_dir/${cohortname} $results_dir
done