results_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
cohort_dir4="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/04/"
#rsync -av /projects/MRC-IEU/research/data/godmc/_devs/GODMC_Analysis/data/sftp/GoDMC/*/*04.tgz $cohort_dir4
cd $cohort_dir4
#gunzip *tgz


rm $results_dir/cohortsizeslambda4.txt
touch $results_dir/cohortsizeslambda4.txt

cd $results_dir

for cohort in ${cohort_dir4}*04.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	echo $cohortname
	mkdir -p $cohort_dir4/${cohortname}
	cd $cohort_dir4/${cohortname}

	#if [ $cohortname == "TwinsUK_04" ]; then
#		samplesize="843"
#		lambda="1.02586328761518"
#	else
    	tar -xvf ../${cohortname}.tar results/04/logs_f/log.txt    
#	

	    #samplesize=`grep founders results/16/logs_b/log.txt |head -n1 |cut -d" " -f5`
     	 samplesize=`grep "people pass filters and QC" results/04/logs_f/log.txt|head -n1 |cut -d" " -f4`
     	 lambda=`grep "lambda value for GWAS:" results/04/logs_f/log.txt|head -n1 |cut -d" " -f5`

#		fi

    echo $cohortname $samplesize $lambda>>$results_dir/cohortsizeslambda4.txt
    
	rm -r $cohort_dir4/${cohortname}
done

samplesize="843"
lambda="1.02586328761518"
cohortname="TwinsUK"

echo $cohortname $samplesize $lambda>>$results_dir/cohortsizeslambda4.txt

mv $results_dir/cohortsizeslambda4.txt /panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda4.txt

cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/01_meta_analysis_16

#extract results files
##results_dir="/projects/MRC-IEU/groups/godmc/meta-analysis/input"
results_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/input"
cohort_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/16/"


rm $results_dir/cohortsizeslambda.txt
touch $results_dir/cohortsizeslambda.txt

cd $results_dir

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	echo $cohortname
	mkdir -p $cohort_dir/${cohortname}
	cd $cohort_dir/${cohortname}

	#if [ $cohortname == "NTR_16" ]; then
#		tar -zxvf ../${cohortname}.tgz results/16/logs_b/log.txt
 #   else
    	tar -xvf ../${cohortname}.tar results/16/logs_b/log.txt    
#	fi

	if [ $cohortname == "BASICMAR_16" ]; then
		samplesize="511"
    else
    	#samplesize=`grep founders results/16/logs_b/log.txt |head -n1 |cut -d" " -f5`
     	 samplesize=`grep "people pass filters and QC" results/16/logs_b/log.txt|tail -n1 |cut -d" " -f4`
     	 lambda=`grep "lambda value for GWAS:" results/16/logs_b/log.txt|tail -n1 |cut -d" " -f5`
	fi

    echo $cohortname $samplesize $lambda>>$results_dir/cohortsizeslambda.txt
    
	rm -r $cohort_dir/${cohortname}
done

mv $results_dir/cohortsizeslambda.txt /panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/cohortsizeslambda.txt