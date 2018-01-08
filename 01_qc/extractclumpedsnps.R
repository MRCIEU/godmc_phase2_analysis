path="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis"

load(paste0(path,"/results/16/16_clumped.rdata"))
dim(clumped)
#[1] 342722     28

clumped<-clumped[clumped$pval<1e-14,]
dim(clumped)
#[1] 288797     28

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
#435391

#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509

nrow(clumped)
#288797
clumped<-clumped[which(clumped$cpg%in%retaincpg),]
nrow(clumped)
#226205

indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt",sep="\t",he=F)
w<-which(clumped$snp%in%indels[,1])
if(length(w)>0){
clumped<-clumped[-w,]
}

write.table(unique(clumped$snp),"clumpedsnps.txt",sep="\t",row.names=F,quote=F,col.names=F)

