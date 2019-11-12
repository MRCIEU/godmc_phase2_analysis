library(data.table)
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)


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

r<-separate_rows(r,genesUniq,sep = ";")
nrow(r) #541872

r2<-fread("EPIC.hg38.manifest.gencode.v22.tsv.gz")
nrow(r2) #865918
m<-match(rm,r2$probeID)
r2<-r2[-na.omit(m),]
nrow(r2) #805311

r2<-separate_rows(r2,genesUniq,sep = ";")
nrow(r2) #1007473

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

probecount<-group_by(r, genesUniq) %>% summarise(nprobes=n(),nmqtls=sum(mqtls,na.rm=T),ncismqtls=sum(cismqtls,na.rm=T),ntransmqtls=sum(transmqtls,na.rm=T))
nrow(probecount)
#33480

probecount2<-group_by(r2, genesUniq) %>% summarise(nprobesEPIC=n())
nrow(probecount2)
#[1] 44633

p<-merge(probecount,probecount2,by.x="genesUniq",by.y="genesUniq",all=T)
nrow(p)
#44709
p[p=="NA"] <- NA
w<-which(is.na(p$nprobes)
sum(p$nmqtls,na.rm=T) #340977
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

w<-which(is.na(p$genesUniq))
pl1<-ggplot(p[-w,], aes(x=bins, y=value)) +
geom_boxplot() +
facet_wrap(~what,ncol=1,scales="free_y", strip.position = "left", labeller = as_labeller(c(nmqtl = "number of mQTL", ratio = "EPIC/450k ratio",n_genes="Number of genes") ) )  +
ylab(NULL)+
theme(axis.text.x=element_text(angle=90, size=10)) +
labs(x="Number of probes per gene")
ggsave(pl1,file="mqtlbygene.pdf",height=10,width=10)


######genecount 450k
genecount<-r %>% select(genesUniq)%>%table%>%table
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
genecount2<-r2 %>% select(genesUniq)%>%table%>%table
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

g<-data.frame(genesUniq="NA",nprobes=g$nprobes,nmqtls=NA,ncismqtls=NA,ntransmqtls=NA,nprobesEPIC=NA,bins=g$bins,ratio=NA,what="ngenes",value=g$ngenes,fill=g$what)

df<-rbind(g,p)
table(df$fill,df$what)

w<-which(df$what!="nmqtl")
df2<-df[w,]

w<-which(is.na(df2$genesUniq))
df2<-df2[-w,]

table(df2$fill,df2$what)

pl1<-ggplot(df2, aes(x=bins, y=value,fill=fill)) +
geom_boxplot() +
facet_wrap(~what,ncol=1,scales="free_y", strip.position = "left", labeller = as_labeller(c(nmqtl_cistrans = "number of mQTL", ratio = "EPIC/450k ratio",ngenes="Number of genes") ) )  +
ylab(NULL)+
theme(axis.text.x=element_text(angle=90, size=10)) +
#scale_fill_discrete(name="array") +
theme(legend.title = element_blank())+
labs(x="Number of probes per gene")
ggsave(pl1,file="mqtlbycoverage.pdf",height=10,width=10)

summary(lm(p1$nmqtls~p1$nprobes))

w<-which(df2$what%in%c("nmqtl_cis","nmqtl_trans"))
pl1<-ggplot(df2[w,], aes(x=bins, y=value,fill=fill)) +
geom_boxplot(outlier.shape= NA) +
facet_wrap(~what,ncol=1,scales="free_y", strip.position = "left", labeller = as_labeller(c(nmqtl_cis = "number of cis mQTL", nmqtl_trans = "number of trans mQTL" ) ))  +
ylab(NULL)+
theme(axis.text.x=element_text(angle=90, size=10)) +
#scale_fill_discrete(name="array") +
theme(legend.title = element_blank())+
labs(x="Number of probes per gene")
ggsave(pl1,file="mqtlbynumberofgenes.pdf",height=10,width=10)

#####
w<-which(df2$what%in%c("nmqtl_cis"))
pl2<-ggplot(df2[w,], aes(x=nprobes)) +
geom_histogram() +
theme(axis.text.x=element_text(angle=90, size=10)) +
labs(x="Number of probes per gene")
ggsave(pl2,file="numberofprobespergene.pdf",height=6,width=10)

w<-which(df2$what%in%c("nmqtl_cis"))
pl2<-ggplot(df2[w,], aes(x=nmqtls)) +
geom_histogram() +
theme(axis.text.x=element_text(angle=90, size=10)) +
labs(x="Number of mqtl per gene")
ggsave(pl2,file="numberofmqtlspergene.pdf",height=6,width=10)

df3<-df2[w,]
df2.plot <- df2 %>%
group_by(bins) %>%
summarize(median.bin = median(nprobes, na.rm=TRUE),median.cis.mqtl=median(ncismqtls, na.rm=TRUE),median.trans.mqtl=median(ntransmqtls, na.rm=TRUE))

summary(lm(df2.plot$median.cis.mqtl~df2.plot$median.bin))

#Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -7.45163    3.98004  -1.872   0.0766 .  
#df2.plot$median.bin  0.75963    0.03127  24.296 9.03e-16 ***

summary(lm(df2.plot$median.trans.mqtl~df2.plot$median.bin))
#Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -0.530955   1.027345  -0.517    0.611    
#df2.plot$median.bin  0.046552   0.008071   5.768 1.47e-05 ***

#number of mQTL from probes with NA genes.
#calculate slope for nmQTL vs number of probes
#mQTL~nprobes

For each probe, you will find an increase of 0.75 mQTL.

summary(lm(p1$nmqtls~p1$nprobes))



#Call:
#lm(formula = p1$nmqtls ~ p1$nprobes)

#Residuals:
#     Min       1Q   Median       3Q      Max 
#-120.814   -3.632    1.153    3.153  193.584 

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -2.751652   0.047016  -58.53   <2e-16 ***
#p1$nprobes   0.799270   0.000147 5437.73   <2e-16 ***

summary(lm(p1$ncismqtls~p1$nprobes))

#Call:
#lm(formula = p1$ncismqtls ~ p1$nprobes)

#Residuals:
#     Min       1Q   Median       3Q      Max 
#-148.883   -3.404    1.245    3.018  207.708 

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -2.7292183  0.0425608  -64.12   <2e-16 ***
#p1$nprobes   0.7422052  0.0001331 5578.07   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.778 on 33478 degrees of freedom
#  (11229 observations deleted due to missingness)
#Multiple R-squared:  0.9989,  Adjusted R-squared:  0.9989 
#F-statistic: 3.111e+07 on 1 and 33478 DF,  p-value: < 2.2e-16

summary(lm(p1$ntransmqtls~p1$nprobes))

#Call:
#lm(formula = p1$ntransmqtls ~ p1$nprobes)

#Residuals:
#    Min      1Q  Median      3Q     Max 
#-42.668  -0.662  -0.206  -0.005  75.859 

#Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)    
#(Intercept) -2.243e-02  1.189e-02   -1.887   0.0592 .  
#p1$nprobes   5.707e-02  3.718e-05 1534.973   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 2.173 on 33478 degrees of freedom
#  (11229 observations deleted due to missingness)
#Multiple R-squared:  0.986, Adjusted R-squared:  0.986 
#F-statistic: 2.356e+06 on 1 and 33478 DF,  p-value: < 2.2e-16


#one additional probe will increase the number of cismQTL with 0.74 (p<2e-16, se=0.0001331) and the number of trans mQTL with 0.057 (p<2e-16,se = 3.718e-05)

