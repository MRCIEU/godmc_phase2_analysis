library(data.table)

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]

r<-fread("Sugden_MethylationReliability_Data_S1.txt",sep="\t",he=T)
r<-as.data.frame(r)
poor<-r[which(r$Reliability<0.4),]
w<-which(poor[,1]%in%retaincpg)
length(w) #299287
nrow(poor) #337550

good<-r[which(r$Reliability>0.75),]
nrow(good) #28249
i<-intersect(good[,1],retaincpg)
length(i) #19199

