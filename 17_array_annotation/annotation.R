library(data.table)
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
ann450 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
length(unique(ann450$Name))
#485512

library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
length(unique(annEPIC$Name))
#866836

ov<-intersect(annEPIC$Name,ann450$Name)
#[1] 453093

w<-which(annEPIC$Name%in%ov)
annEPIC2<-annEPIC[-w,]
head(data.frame(annEPIC2))

w<-which(annEPIC2$Islands_Name%in%ann450$Islands_Name)

isl<-unique(ann450$Islands_Name)
isl<-isl[which(isl%in%""==F)]

#non 450k CpGs on 450k island
annEPIC3<-annEPIC2[which(annEPIC2$Islands_Name%in%isl),]
#90103

rf<-unique(ann450$Regulatory_Feature_Name)
rf<-rf[which(rf%in%""==F)]

#non 450k CpGs on 450 regulatory feature
annEPIC4<-annEPIC2[which(annEPIC2$Regulatory_Feature_Name%in%rf),]
#35513

ov<-intersect(annEPIC3$Name,annEPIC4$Name)
#29187

tepic<-annEPIC[annEPIC$Name%in%ov,]

df1<-data.frame(annEPIC[which(annEPIC$Regulatory_Feature_Name%in%c("15:74286538-74287451")),c("chr","pos","Name")])
df2<-data.frame(ann450[which(ann450$Regulatory_Feature_Name%in%c("15:74286538-74287451")),c("chr","pos","Name")])

intersect(df1$Name, df2$Name)

length(unique(tepic$Regulatory_Feature_Name))
#[1] 12423
length(unique(tepic$Islands_Name))
#11194

positions = data.frame (epic=annEPIC2$Name,chr=annEPIC2$chr,position=annEPIC2$pos)
intervals = data.frame (a450k=ann450$Name,chr=ann450$chr,start=ann450$pos,end=ann450$pos)

genomic.nearest <- function(positions, intervals) {
    events <- rbind(data.frame(chr=intervals$chr,
                               position=intervals$start,
                               type="start",
                               id=1:nrow(intervals)),
                    data.frame(chr=intervals$chr,
                               position=intervals$end,
                               type="end",
                               id=1:nrow(intervals)),
                    data.frame(chr=positions$chr,
                               position=positions$position,
                               type="position",
                               id=1:nrow(positions)))
    events <- events[order(events$chr, events$position, decreasing=F),]

    before <- (1:nrow(events))-1
    before[1] <- NA
    for (position.idx in which(events$type == "position" & events$type[before] == "position"))
        if (position.idx > 1)
            before[position.idx] <- before[position.idx-1] ## 1s                                                                                                                                            

    after <- (1:nrow(events))+1
    after[length(after)] <- NA
    for(position.idx in rev(which(events$type == "position" & events$type[after] == "position")))
        if (position.idx < nrow(events))
            after[position.idx] <- after[position.idx+1]


    dist.before <- events$position - events$position[before]
    dist.before[which(events$chr != events$chr[before])] <- NA
    dist.after <- events$position[after] - events$position
    dist.after[which(events$chr != events$chr[after])] <- NA
    events$nearest <- NA
    idx <- which(!is.na(dist.before) & (is.na(dist.after) | dist.before <= dist.after))
    events$nearest[idx] <- before[idx]
    idx <- which(!is.na(dist.after) & (is.na(dist.before) | dist.before > dist.after))
    events$nearest[idx] <- after[idx]

    events$nearest <- events$id[events$nearest]
    positions$nearest <- NA
    idx <- which(events$type == "position")
    positions$nearest[events$id[idx]] <- events$nearest[idx]
    positions$nearest
}

idx<-genomic.nearest(positions,intervals)
length(idx)
#[1] 413743

df<-data.frame(positions,intervals[idx,])
df$dist<-abs(df$position-df$start)
hist(df$dist)

quantile(df$dist)
#      0%      25%      50%      75%     100% 
#       2      459     2372     8609 30463967 

length(which(df$dist<10000)) #78%
#321373

length(which(df$dist<1000))/nrow(df)
#148246
#[1] 0.3583046

length(which(df$dist<100))/nrow(df)
#[1] 0.1152092

#78% (n=321373) of non450k probes are within 10k of 450k probe
#36% (n=148246) of non450 probes are within 1k of 450k probe
#12% (n=47667) of non 450 probes are within 100bp of 450k probe
#~12,000 regulatory features contain non450k probes that map to same regulatory feature as 450k probe 

####################################################################################################

#x-axis cpg density
#y-axis number of mQTL
load("~/repo/godmc_phase2_analysis/07_enrichments/mean_allcpgs.Robj")

r<-fread("hm450.hg38.manifest.gencode.v22.tsv.gz")
nrow(r) #485577
rm<-which(r$probeID%in%df.all$cpg==F)
rm<-r$probeID[rm]
m<-match(rm,r$probeID)
r<-r[-m,]
nrow(r) #420509


which(is.na(r$genesUniq))

r<-separate_rows(r,geneNames,transcriptTypes,sep = ";")
nrow(r) #2068927

r<-unique(r)
nrow(r)
#542090

r.rm<-r[which(is.na(r$transcriptTypes)),]

genes<-unlist(c("protein_coding",unique(r2[grep("_gene",as.character(r2$transcriptTypes),fixed=T),"transcriptTypes"])))
r<-r[which(r$transcriptTypes%in%genes),]
length(unique(r$geneNames)) #18993
nrow(r) #386057
length(unique(r$probeID)) #331884
nrow(r) #386057


r2<-fread("EPIC.hg38.manifest.gencode.v22.tsv.gz")
nrow(r2) #865918
m<-match(rm,r2$probeID)
r2<-r2[-na.omit(m),]
nrow(r2) #805311

r2<-separate_rows(r2,geneNames,transcriptTypes,sep = ";")
nrow(r2) #3735330

r2<-unique(r2)
nrow(r2)
#1007862

r.rm2<-r2[which(is.na(r2$transcriptTypes)),]

genes<-unlist(c("protein_coding",unique(r2[grep("_gene",as.character(r2$transcriptTypes),fixed=T),"transcriptTypes"])))
r2<-r2[which(r2$transcriptTypes%in%genes),]
length(unique(r2$geneNames)) #19543
nrow(r2) #659907


load("~/repo/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped<-clumped[which(clumped$cis==TRUE & clumped$pval<1e-8 | clumped$cis==FALSE & clumped$pval<1e-14),]
mqtlcount <- group_by(clumped, cpg) %>% summarise(mqtls=n())
clumped.cis<-clumped[which(clumped$cis==TRUE & clumped$pval<1e-8),]
clumped.trans<-clumped[which(clumped$cis==FALSE & clumped$pval<1e-14),]

ciscount<-group_by(clumped.cis, cpg) %>% summarise(cismqtls=n())
transcount<-group_by(clumped.trans, cpg) %>% summarise(transmqtls=n())

r<-merge(r,mqtlcount,by.x="probeID",by.y="cpg",all=T)
r<-merge(r,ciscount,by.x="probeID",by.y="cpg",all=T)
r<-merge(r,transcount,by.x="probeID",by.y="cpg",all=T)

r.rm<-merge(r.rm,mqtlcount,by.x="probeID",by.y="cpg",all.x=T)
r.rm<-merge(r.rm,ciscount,by.x="probeID",by.y="cpg",all.x=T)
r.rm<-merge(r.rm,transcount,by.x="probeID",by.y="cpg",all.x=T)


probecount<-group_by(r, geneNames) %>% summarise(nprobes=n(),nmqtls=sum(mqtls,na.rm=T),ncismqtls=sum(cismqtls,na.rm=T),ntransmqtls=sum(transmqtls,na.rm=T))
nrow(probecount)
#18994

probecount.rm<-group_by(r.rm, geneNames) %>% summarise(nprobes=n(),nmqtls=sum(mqtls,na.rm=T),ncismqtls=sum(cismqtls,na.rm=T),ntransmqtls=sum(transmqtls,na.rm=T))
nrow(probecount.rm)
# A tibble: 1 x 5
#  geneNames nprobes nmqtls ncismqtls ntransmqtls
#  <chr>       <int>  <int>     <int>       <int>
#1 NA          58356  46654     43324        3330

probecount2<-group_by(r2, geneNames) %>% summarise(nprobesEPIC=n())
nrow(probecount2)
#[1] 19543

probecount.rm2<-group_by(r.rm2, geneNames) %>% summarise(nprobesEPIC=n())
p.rm<-merge(probecount.rm,probecount.rm2,by.x="geneNames",by.y="geneNames",all=T)

p<-merge(probecount,probecount2,by.x="geneNames",by.y="geneNames",all=T)
nrow(p)
#19548
p[p=="NA"] <- NA
w<-which(is.na(p$nprobes)) #554
sum(p$nmqtls,na.rm=T) #300831
#p[which(is.na(p$nmqtls)),"nmqtls"]<-0
#p[which(is.na(p$ncismqtls)),"ncismqtls"]<-0
#p[which(is.na(p$ntransmqtls)),"ntransmqtls"]<-0
p[which(is.na(p$nprobes)),"nprobes"]<-0
p[which(is.na(p$nprobesEPIC)),"nprobesEPIC"]<-0

#breaks <- c(seq(0,200,10),seq(200,1100,100))
breaks<-c(0,seq(1,200,10),201,max(p$nprobes))
#breaks[2]<-1
bins <- cut(p$nprobes, breaks, include.lowest = T, right=FALSE)
bins<-do.call("rbind",strsplit(as.character(bins),split=","))
bins<-gsub("[","",bins,fixed=T)
bins<-gsub("]","",bins,fixed=T)
bins<-gsub(")","",bins)
bins[,2]<-as.numeric(bins[,2])-1

p<-data.frame(p[,1:6],bins=as.character(paste0(bins[,1],"-",bins[,2])))
p$bins<-as.character(p$bins)


w<-which(p$nprobes>200)
p$bins<-as.character(p$bins)
p$bins[w]<-">200"
w<-which(p$nprobes=="0")
p$bins[w]<-"0"

o<-order(as.numeric(unique(bins[,1])))
p$bins <- factor(p$bins, levels = c(unique(p$bins)[o]))

p$ratio<-p$nprobesEPIC/p$nprobes
p[which(is.na(p$ratio)),]
w<-which(p$ratio=="Inf")
p$ratio[w]<-p$nprobesEPIC[w]

p1<-p
p1$what<-"nmqtl"
p1$value<-p$nmqtls
p1$fill<-"450k"

p2<-p
p2$what<-"ratio"
p2$value<-p$ratio
p2$fill<-"450k"

p3a<-p
p3a$what<-"nmqtl_cis"
p3a$value<-p$ncismqtls
p3a$fill<-"cis"

p3b<-p
p3b$what<-"nmqtl_trans"
p3b$value<-p$ntransmqtls
p3b$fill<-"trans"

p<-rbind(p1,p2,p3a,p3b)

p0<-p[which(p$nprobes==0&p$what=="ratio"),]
data.frame(table(p0$nprobesEPIC))

p[which(p$what=="nmqtl"&is.na(p$geneNames)),]
#      geneNames nprobes nmqtls ncismqtls ntransmqtls nprobesEPIC bins ratio
#19548      <NA>   49319  71635     65873        5762           0 >200     0
#       what value fill
#19548 nmqtl 71635 450k

w<-which(is.na(p$geneNames))
pl1<-ggplot(p[-w,], aes(x=bins, y=value)) +
geom_boxplot(outlier.shape= NA) +
facet_wrap(~what,ncol=1,scales="free_y", strip.position = "left", labeller = as_labeller(c(nmqtl = "number of mQTL", ratio = "EPIC/450k ratio") ) )  +
ylab(NULL)+
theme(axis.text.x=element_text(angle=90, size=10)) +
labs(x="Number of probes per gene")
ggsave(pl1,file="mqtlbygene.pdf",height=10,width=10)


######genecount 450k
genecount<-r %>% select(geneNames)%>%table%>%table
nrow(genecount)
genecount<-data.frame(genecount)
names(genecount)<-c("nprobes","ngenes")

g0<-data.frame(nprobes=0,ngenes=NA)
genecount<-rbind(g0,genecount)


genecount$nprobes<-as.numeric(as.character(genecount$nprobes))
breaks<-c(0,seq(1,200,10),201,max(genecount$nprobes))
#breaks[2]<-1
bins <- cut(genecount$nprobes, breaks, include.lowest = T, right=FALSE)
bins<-do.call("rbind",strsplit(as.character(bins),split=","))
bins<-gsub("[","",bins,fixed=T)
bins<-gsub("]","",bins,fixed=T)
bins<-gsub(")","",bins)
bins[,2]<-as.numeric(bins[,2])-1

g<-data.frame(genecount,bins=as.character(paste0(bins[,1],"-",bins[,2])))
g$bins<-as.character(g$bins)
w<-which(g$nprobes>200)
g$bins<-as.character(g$bins)
g$bins[w]<-">200"
w<-which(g$nprobes==0)
g$bins[w]<-"0"
g$what<-"450k"

o<-order(as.numeric(unique(bins[,1])))
g$bins <- factor(g$bins, levels = c(unique(g$bins)[o]))
levels(g$bins)
#####genecount EPIC
genecount2<-r2 %>% select(geneNames)%>%table%>%table
nrow(genecount2)

genecount2<-data.frame(genecount2)
names(genecount2)<-c("nprobes","ngenes")

g0<-data.frame(nprobes=0,ngenes=NA)
genecount2<-rbind(g0,genecount2)

genecount2$nprobes<-as.numeric(as.character(genecount2$nprobes))
breaks<-c(0,seq(1,200,10),201,max(genecount2$nprobes))
#breaks[2]<-1
bins <- cut(genecount2$nprobes, breaks, include.lowest = T, right=FALSE)
bins<-do.call("rbind",strsplit(as.character(bins),split=","))
bins<-gsub("[","",bins,fixed=T)
bins<-gsub("]","",bins,fixed=T)
bins<-gsub(")","",bins)
bins[,2]<-as.numeric(bins[,2])-1

g2<-data.frame(genecount2,bins=as.character(paste0(bins[,1],"-",bins[,2])))
g2$bins<-as.character(g2$bins)
w<-which(g2$nprobes>200)
g2$bins<-as.character(g2$bins)
g2$bins[w]<-">200"
w<-which(g2$nprobes==0)
g2$bins[w]<-"0"

g2$what<-"EPIC"

o<-order(as.numeric(unique(bins[,1])))
g2$bins <- factor(g2$bins, levels = c(unique(as.character(g2$bins))[o]))
levels(g2$bins)

###
g<-rbind(g,g2)
#g[which(g$bins%in%c("0","1-10",">200")),]

##combine probe count and gene count for plotting
#g<-data.frame(genesUniq="NA",nprobes=g$nprobes,nmqtls=NA,ncismqtls=NA,ntransmqtls=NA,nprobesEPIC=NA,bins=g$bins,ratio=NA,what="ngenes",value=g$ngenes,fill=g$what)

#df<-rbind(g,p)
#table(df$fill,df$what)

#w<-which(df$what!="nmqtl")
#df2<-df[w,]

#w<-which(is.na(df2$genesUniq))
#df2<-df2[-w,]

#table(df2$fill,df2$what)

pl1<-ggplot(g, aes(x=bins, y=ngenes,fill=what)) +
geom_boxplot() +
#facet_wrap(~what,ncol=1,scales="free_y")  +
ylab(NULL)+
theme(axis.text.x=element_text(angle=90, size=10)) +
#scale_fill_discrete(name="array") +
theme(legend.title = element_blank())+
labs(x="Number of probes per gene")
ggsave(pl1,file="probespergene.pdf",height=6,width=10)

summary(lm(p1$nmqtls~p1$nprobes))

w<-which(p$what%in%c("nmqtl_cis","nmqtl_trans")&!is.na(p$geneNames))
pl1<-ggplot(p[w,], aes(x=bins, y=value,fill=fill)) +
geom_boxplot(outlier.shape= NA) +
facet_wrap(~what,ncol=1,scales="free_y", strip.position = "left", labeller = as_labeller(c(nmqtl_cis = "number of cis mQTL", nmqtl_trans = "number of trans mQTL" ) ))  +
ylab(NULL)+
theme(axis.text.x=element_text(angle=90, size=10)) +
#scale_fill_discrete(name="array") +
theme(legend.title = element_blank())+
labs(x="Number of probes per gene")
ggsave(pl1,file="mqtlbynumberofgenes.pdf",height=10,width=10)

#####
w<-which(p$what%in%c("nmqtl_cis")&!is.na(p$geneNames))
pl2<-ggplot(p[w,], aes(x=nprobes)) +
geom_histogram() +
theme(axis.text.x=element_text(angle=90, size=10)) +
labs(x="Number of probes per gene")
ggsave(pl2,file="numberofprobespergene.pdf",height=6,width=10)

w<-which(p$what%in%c("nmqtl_cis")&!is.na(p$geneNames))
pl2<-ggplot(p[w,], aes(x=nmqtls)) +
geom_histogram() +
theme(axis.text.x=element_text(angle=90, size=10)) +
labs(x="Number of mqtl per gene")
ggsave(pl2,file="numberofmqtlspergene.pdf",height=6,width=10)

table(p$what)

#      nmqtl   nmqtl_cis nmqtl_trans       ratio 
#      19548       19548       19548       19548 

w<-which(p$what%in%c("nmqtl")&!is.na(p$geneNames)&p$bins!="0")

df3<-p[w,]
df2.plot <- df3 %>%
group_by(bins) %>%
summarize(median.bin = median(nprobes, na.rm=TRUE),median.cis.mqtl=median(ncismqtls, na.rm=TRUE),median.trans.mqtl=median(ntransmqtls, na.rm=TRUE))

summary(lm(df2.plot$median.cis.mqtl~df2.plot$median.bin))

#Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -7.06326    3.87338  -1.824    0.084 .  
#df2.plot$median.bin  0.74133    0.03051  24.300    9e-16 ***

summary(lm(df2.plot$median.trans.mqtl~df2.plot$median.bin))
#Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -0.29548    0.28823  -1.025    0.318    
#df2.plot$median.bin  0.03930    0.00227  17.310 4.32e-13 ***

#number of mQTL from probes with NA genes.
#calculate slope for nmQTL vs number of probes
#mQTL~nprobes

#For each probe, you will find an increase of 0.74 mQTL.

head(p[is.na(p$geneNames),])
#      geneNames nprobes nmqtls ncismqtls ntransmqtls nprobesEPIC bins ratio
#19548      <NA>   49319  71635     65873        5762           0 >200     0
#39096      <NA>   49319  71635     65873        5762           0 >200     0
#58644      <NA>   49319  71635     65873        5762           0 >200     0
#78192      <NA>   49319  71635     65873        5762           0 >200     0
#             what value  fill
#19548       nmqtl 71635  450k
#39096       ratio     0  450k
#58644   nmqtl_cis 65873   cis
#78192 nmqtl_trans  5762 trans





