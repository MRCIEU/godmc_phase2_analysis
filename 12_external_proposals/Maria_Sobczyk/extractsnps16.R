res.out<-data.frame()
res.out2<-data.frame()

#rsid
snps<-read.csv("/panfs/panasas01/shared-godmc/database_files/snps.csv")
r<-read.table("paternoster_2015_index_snps_sorted_3Mbp.snps")
w<-which(snps$rsid%in%r$V1)
snps<-snps[w,]

#chr:pos
r2<-read.table("paternoster_2015_index_snps_sorted_3Mbp.godmc",he=T)


for (i in 1:962){
cat(i,"\n")
p<-paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_",i,".rdata")
load(p)
res1<-res[which(res$snp%in%snps[,1]),]
res.out<-rbind(res.out,res1)

res2<-res[which(res$snp%in%r2$name),]
res.out2<-rbind(res.out2,res2)

}
write.table(res.out,"~/repo/godmc_phase2_analysis/12_external_proposals/Maria_Sobczyk/snps_36cohorts16_rsid.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(res.out2,"~/repo/godmc_phase2_analysis/12_external_proposals/Maria_Sobczyk/snps_36cohorts16_chrpos.txt",sep="\t",col.names=T,row.names=F,quote=F)


i<-intersect(unique(res.out$snp),snps$name)
length(i)
length(unique(res.out$snp))

i<-intersect(unique(res.out2$snp),r2$V1)
length(i)