#extract results files
results_dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/16"
outdir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/control"
cohort_dir="../data/16/"
#metal_in="16_${i}.in"

cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/01_meta_analysis_16

for cohort in ${cohort_dir}*16.tar
do

	echo $cohort
	cohortname=$(basename "$cohort" .tar)
	cohortname2=$(basename "$cohort" _16.tar)
	echo $cohortname
	mkdir -p $results_dir/${cohortname}
	
	mkdir -p $outdir/${cohortname}
	cd $results_dir/${cohortname}
	tar -xvf ../${cohortname}.tar results/16/control/cg07959070.assoc.linear.gz
	echo "CHR" "MARKERNAME" "POS" "EA" "MODEL" "N" "STAT" "BETA" "P" |gzip -c > header.gz
	zcat header.gz $results_dir/${cohortname}/results/16/control/cg07959070.assoc.linear.gz | gzip -c  >$outdir/${cohortname}.control.gz
	zcat $outdir/${cohortname}.control.gz | sed 's/ \+/ /g' |sed 's/^ //g' | perl -pe 's/ /\t/g'|gzip -c > $outdir/${cohortname}.control.gz2
	mv $outdir/${cohortname}.control.gz2 $outdir/${cohortname}.control.gz
	cd ..
	rm -r $results_dir/${cohortname}
	cd $results_dir
	perl ~/repo/godmc_phase1_analysis/02.extract_counts/scripts/join_file.pl -i "$outdir/${cohortname}.control.gz,TAB,1 /panfs/panasas01/shared-godmc/godmc_phase2_analysis/scratch/frq/${cohortname2}.frq.tmp,TAB,1" -o $outdir/${cohortname}.control.merged.tmp -a 1
	awk '{ if(NR>1 && $7!="NA") { print $0, $7/$8 } else {print "CHR","MARKERNAME","POS","A1","TEST","NMISS","BETA","TSTAT","P","CHR","SNP","EA","NEA","EAF","N","SE"}}' < $outdir/${cohortname}.control.merged.tmp >$outdir/${cohortname}.control.merged
	
    
	#rm $outdir/${cohortname}.control.merged.tmp

done

##
#se=beta/TSTAT

 #CHR                     SNP         BP   A1       TEST    NMISS       BETA         STAT            P 
 #  1         chr1:768448:SNP     768448    A        ADD     1774    0.04856        1.183       0.2368
 #  1         chr1:768448:SNP     768448    A       COV1     1774     0.9549        36.23   1.447e-215
 #  1         chr1:769963:SNP     769963    A        ADD     1774    0.04856        1.183       0.2368
 #  1         chr1:769963:SNP     769963    A       COV1     1774     0.9549        36.23   1.447e-215
 #  1         chr1:985797:SNP     985797    G        ADD     1774   -0.03614      -0.4623       0.6439
 #  1         chr1:985797:SNP     985797    G       COV1     1774     0.9549        36.22   1.977e-215
 #  1         chr1:993360:SNP     993360    C        ADD     1763     0.0422        1.632        0.103
 #  1         chr1:993360:SNP     993360    C       COV1     1763     0.9542        36.21   4.967e-215
 #  1         chr1:995481:SNP     995481    T        ADD     1767    0.04013        1.547        0.122