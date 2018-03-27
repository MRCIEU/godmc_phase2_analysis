library(data.table)
library(ggplot2)
library(dplyr)

#ov.all<-data.frame()
#for (i in 1:23){
#cat(i,"\n")
#load(paste0("chr",i,".Robj"))
#ov.all<-rbind(ov.all,ov)
#}
#save(ov.all,file="tfbsbycpg.Robj")

load("tfbsbycpg.Robj")
path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/encode_tfbs"
l<-list.files(path)
l2<-gsub("tfbs_","",l)
names(ov.all)[6:ncol(ov.all)]<-l2
save(ov.all,file="tfbsbycpg.Robj")

nvars<-names(ov.all)
#ov.all2<-ov.all[,seq(2,6:ncol(ov.all)),with=FALSE]
#row.names(ov.all2)<-ov.all[,1]
#out<-ov.all%>%group_by(paste0(nvars))%>%summarise(mean=mean(meancpg))
#out<-ov.all[, lapply(.SD,mean(meancpg)), by=ov.all[,6]]
#out<-ov.all[, lapply(.SD,mean), by=eval(colnames(ov.all)[6:ncol(ov.all)])]
#.SD is the (S)ubset of (D)ata excluding group columns
ov.all<-data.frame(ov.all)
df.out <- data.frame()
    for (i in 6:ncol(ov.all)) {
        cat(i,"\n")
        results <- ov.all %>% group_by(ov.all[,i]) %>% summarize(mean=mean(meancpg))
        df<-data.frame(tfbs=results[,1],mean=results[,2],antibody=nvars[i])
        df.out<-rbind(df.out,df)
    }

save(df.out,file="tfbsbycpgmean.Robj")

df.out <- data.frame()
    for (i in 6:ncol(ov.all)) {
        cat(i,"\n")
        results <- ov.all %>% group_by(ov.all[,i]) %>% summarize(mean=mean(sdcpg))
        df<-data.frame(tfbs=results[,1],sd=results[,2],antibody=nvars[i])
        df.out<-rbind(df.out,df)
    }

save(df.out,file="tfbsbycpgsd.Robj")