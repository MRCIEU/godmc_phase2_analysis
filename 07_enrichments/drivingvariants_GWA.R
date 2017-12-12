#awk '{print "chr" $5}' <garfield.test.ambivalentmqtl_sds.out.significant.annotations.1e-14.0.010050256375866.variants|sed 's/[^(]*$//'|perl -pe 's/\(//g' >snps.txt
#awk '{print "chr" $5}' <garfield.test.transmqtl_sds.out.significant.annotations.1e-14.0.0100477631605532.variants|sed 's/[^(]*$//'|perl -pe 's/\(//g' >snps.txt
#ls /mnt/storage/private/mrcieu/research/GODMC_Analysis/godmc_phase2_analysis/data/gwas

cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments
cat /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/transmqtl_sds/snps.txt /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/ambivalentmqtl_sds/snps.txt |sort -u >transsdssnps.txt
awk '{print "chr"$1":SNP"}' <transsdssnps.txt >transsdssnps.plink.txt
plink --bfile /panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered --extract /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments/transsdssnps.plink.txt --make-bed --out test --maf 0.01
plink --bfile test --indep-pairwise 1000kb 5 0.1 --out indep --maf 0.01

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