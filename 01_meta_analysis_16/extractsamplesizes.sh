cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/01_meta_analysis_16

#extract results files
##results_dir="/projects/MRC-IEU/groups/godmc/meta-analysis/input"
results_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
cohort_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/16/"
rm $results_dir/cohortsizes.txt
touch $results_dir/cohortsizes.txt
#metal_in="16_${i}.in"
cd $results_dir

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	echo $cohortname
	mkdir -p $cohort_dir/${cohortname}
	cd $cohort_dir/${cohortname}

	if [ $cohortname == "NTR_16" ]; then
		tar -zxvf ../${cohortname}.tgz results/16/logs_b/log.txt
    else
    	tar -xvf ../${cohortname}.tar results/16/logs_b/log.txt    
	fi

	if [ $cohortname == "BASICMAR_16" ]; then
		samplesize="529"
    else
    	#samplesize=`grep founders results/16/logs_b/log.txt |head -n1 |cut -d" " -f5`
     	 samplesize=`grep "people pass filters and QC" results/16/logs_b/log.txt|tail -n1 |cut -d" " -f4`
	fi


    echo $cohortname $samplesize >>$results_dir/cohortsizes.txt 
	rm -r $cohort_dir/${cohortname}
done
