# script to reformat GWAS p values

#The files must contain no header and have genomic
#position (build 37) in the first column and P-value from association analysis in the second column
#using a space as a column delimiter. 

PVALDIR=/mnt/data1/programs/GARFIELD/garfield-data/pval

mkdir -p $PVALDIR/ALS_LMM
mkdir -p $PVALDIR/ALS_Meta

cd /mnt/data1/reference_files/ALS/
for chr in {1..22};
	do cut -d" " -f 1,3 als.sumstats.lmm.chr${chr}.txt | tail -n +2 | sort > $PVALDIR/ALS_LMM/chr${chr};
done

for chr in {1..22};
	do cut -d" " -f 1,3 als.sumstats.meta.chr${chr}.txt | tail -n +2 | sort > $PVALDIR/ALS_Meta/chr${chr};
done

mkdir -p $PVALDIR/Alzheimers

cd /mnt/data1/reference_files/Alzheimers

awk 'BEGIN { for (chr = 1; chr <= 22; ++i){ if($1 == chr); print $0 }}' IGAP_stage_1.txt 

for chr in {1..22};
	do awk '{if ($1 == chr) print $2,$8}' IGAP_stage_1.txt > $PVALDIR/Alzheimers/chr${chr};
done
	
