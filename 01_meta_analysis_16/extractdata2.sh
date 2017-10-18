#extract results files
#/panfs/panasas01/shared-godmc/sftp/GoDMC/ecarnero/TwinsUK/results/01
cohort_dir="/panfs/panasas01/shared-godmc/results/01/"
results_dir="/projects/MRC-IEU/research/data/godmc/_devs/GODMC_Analysis/data/sftp/GoDMC/"


cd $results_dir
rm /panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortrelatedness.txt
for cohort in ${cohort_dir}*01.tgz
do
echo $cohort


	cohortname=$(basename "$cohort" .tgz)
	echo $cohortname
	mkdir -p $cohort_dir/${cohortname}
	cd $cohort_dir/${cohortname}


	tar -xvf ../${cohortname}.tgz results/01/cohort_descriptives.RData
    tar -xvf ../${cohortname}.tgz config
    related=`grep "related=" config |perl -pe 's/related=//g'|perl -pe 's/"//g'`
    echo $cohortname $related >>/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortrelatedness.txt 
    tar -xvf ../${cohortname}.tgz results/01/methylation_summary.RData
done