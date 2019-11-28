
#!bin/bash

results_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
#don't remove last slash
cohort_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/16/"
outdir="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/01_qc"

mkdir -p $outdir/chr1.assoc
cd $outdir

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	echo $cohortname
	mkdir -p $cohort_dir/${cohortname}

		cd $outdir
		file="chr1_chunks.txt"
		while read -r line
		do
			chunk="$line"
			echo "$chunk"

			cd $cohort_dir/${cohortname}
			tar -xvf ../${cohortname}.tar results/16/results_${chunk}.gz
			zcat results/16/results_${chunk}.gz | fgrep -w -f $outdir/chr1_ids.txt >$outdir/chr1.assoc/${cohortname}.chr1.assoc_${chunk}.txt
		done < $file

###for i in `seq 1 25`; 
###do
###	tar -xvf ../${cohortname}.tar results/16/results_${i}.gz
###done
rm -r $results_dir/${cohortname}
mv $cohort_dir/${cohortname} $results_dir
done
