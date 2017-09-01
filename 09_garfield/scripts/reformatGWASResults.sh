# script to reformat GWAS p values

#The files must contain no header and have genomic
#position (build 37) in the first column and P-value from association analysis in the second column
#using a space as a column delimiter. 

PVALDIR=/mnt/data1/programs/GARFIELD/garfield-data/pval

mkdir -p $PVALDIR/ALS_LMM
mkdir -p $PVALDIR/ALS_Meta

cd /mnt/data1/reference_files/ALS/
for chr in {1..22};
	do cut -d" " -f 3,9 --output-delimiter=' ' als.sumstats.lmm.chr${chr}.txt | tail -n +2 | sort -n > $PVALDIR/ALS_LMM/chr${chr};
done

for chr in {1..22};
	do cut -d" " -f 3,9 als.sumstats.meta.chr${chr}.txt | tail -n +2 | sort -n > $PVALDIR/ALS_Meta/chr${chr};
done

mkdir -p $PVALDIR/Alzheimers

cd /mnt/data1/reference_files/Alzheimers

for chr in {1..22};
	do awk '{if ($1 == "'$chr'") print $2,$8}' IGAP_stage_1.txt | sort -n > $PVALDIR/Alzheimers/chr${chr};
done

mkdir -p $PVALDIR/Anorexia

cd /mnt/data1/reference_files/Anorexia

for chr in {1..22};
	do awk '{if ($1 == "'$chr'") print $5,$9}' pgc.an.snp.all.13May2016.txt | sort -n > $PVALDIR/Anorexia/chr${chr};
done

mkdir -p $PVALDIR/Asthma

cd /mnt/data1/reference_files/Asthma/GWAS/
for chr in {1..22};
	do awk '{if ($1 == "'$chr'") print $3,$11}' gabriel_asthma_meta-analysis_36studies_format_repository_NEJM.txt | sort -n > $PVALDIR/Asthma/chr${chr};
done

mkdir -p $PVALDIR/BipolarDisorder

cd /mnt/data1/reference_files/BipolarDisorder/NIMHGWAS/
for chr in {1..22};
	do awk '{if ($1 == "'$chr'") print $2,$10}' Bipolar_disorder_PCAadj_GC_Neff_meta.tbl | sort -n > $PVALDIR/BipolarDisorder/chr${chr};
done

mkdir -p $PVALDIR/CAD_2014
mkdir -p $PVALDIR/MI_2015


cd /mnt/data1/reference_files/Cardiogram
for chr in {1..22};
	do awk '{if ($2 == "'$chr'") print $3,$11}' cad.add.160614.website.txt | sort -n > $PVALDIR/CAD_2014/chr${chr};
done

cd /mnt/data1/reference_files/Cardiogram
for chr in {1..22};
	do awk '{if ($2 == "'$chr'") print $3,$11}' mi.add.030315.website.txt | sort -n > $PVALDIR/MI_2015/chr${chr};
done


