ss<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/descriptives/descriptives.phase2.txt",sep="\t",he=T)
ss<-ss[,c("study","nsamples04")]
ss1<-data.frame(study="DunedinAge38",nsamples04="757")
ss<-rbind(ss,ss1)
cohorts<-read.table("/panfs/panasas01/shared-godmc/scripts/cohorts_noRAINE.txt",sep=" ",header=F)
ss<-ss[which(ss$study%in%cohorts$V2),]
save(ss, file="studies.rdata")

