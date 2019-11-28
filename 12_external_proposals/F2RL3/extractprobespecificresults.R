res.out<-data.frame()
for (i in 1:962){
cat(i,"\n")
p<-paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_cleaned_",i,".rdata")
load(p)
res<-res[which(res$cpg%in%c("cg03636183")),]
res.out<-rbind(res.out,res)
}
save(res.out,file="cg03636183_n36.RData")
