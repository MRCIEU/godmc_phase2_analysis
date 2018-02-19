#awk '{print "chr" $5}' <garfield.test.ambivalentmqtl_sds.out.significant.annotations.1e-14.0.010047802254319.variants|sed 's/[^(]*$//'|perl -pe 's/\(//g' >snps.txt
#awk '{print "chr" $5}' <garfield.test.cismqtl_sds.out.significant.annotations.1e-8.0.0100524490036318.variants|sed 's/[^(]*$//'|perl -pe 's/\(//g' >snps.txt
#ls /mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/data/gwas

cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/13_selection
cat /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/cismqtl_sds/snps.txt /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/ambivalentmqtl_sds/snps.txt |sort -u >cissdssnps.txt
awk '{print $1":SNP"}' <cissdssnps.txt >cissdssnps.plink.txt
plink --bfile /panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered --extract /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/13_selection/cissdssnps.plink.txt --make-bed --out test --maf 0.01
plink --bfile test --indep-pairwise 1000kb 5 0.001 --out indep --maf 0.01

#
library(biomaRt)
listEnsembl()

snps<-read.table("~/repo/godmc_phase2_analysis/07_enrichments/snps.txt")
snps[,1]<-gsub("chr","",snps[,1])
#spl<-strsplit(as.character(snps[,1]),split=":")
#spl<-do.call("rbind",spl)

variation = useEnsembl(biomart="snp", dataset="hsapiens_snp",GRCh=37)
rsid <- getBM(attributes=c('refsnp_id','refsnp_source','chr_name','chrom_start','chrom_end'), filters = 'chromosomal_region', values =snps[,1], mart = variation)

p<-paste(rsid$chr_name,rsid$chrom_end,sep=":")
w<-which(p%in%snps[,1])

p<-paste("chr",spl[,1],":",spl[,2],"-",spl[,2])
