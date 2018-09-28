library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(dplyr)

load("../results/graph.rdata")
df<-data.frame(table(mem$cluster))
clusters<-df[which(df$Freq>=10),1]

load("../results/core_communities_cpg_tophits.rdata")
max(core_communities_cpg_tophits$fdr)
#1
max(core_communities_cpg_tophits$fdr2)
#0.04996989

length(unique(core_communities_cpg_tophits$antibody)) #120
core_communities_cpg_tophits$antibody2<-sub("_[^_]+$", "", core_communities_cpg_tophits$antibody)
length(unique(core_communities_cpg_tophits$antibody2)) #115
core_communities_cpg_tophits<-core_communities_cpg_tophits[which(core_communities_cpg_tophits$userSet%in%clusters),]
length(unique(core_communities_cpg_tophits$antibody2)) #98

pval_lim<-0.001

core_communities_cpg_tophits$Pvalue<-10^(-core_communities_cpg_tophits$pValueLog)
 
    m<-max(core_communities_cpg_tophits$logOddsRatio)+1
    p1<-ggplot(core_communities_cpg_tophits,aes(x=antibody2,y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=userSet))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
    #guides(fill = guide_legend(ncol=20))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Transcription factor", fill="community")
      
    ggsave(p1,file=paste0("../images/epigenetic_TFBS_OR.pdf"),height=10,width=16)

core_communities_cpg_tophits$userSet<-factor(core_communities_cpg_tophits$userSet,levels=clusters)
  p1<-ggplot(core_communities_cpg_tophits,aes(x=factor(userSet),y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=antibody2))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
      guides(fill = guide_legend(nrow=2))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Transcription factor",fill="TFBS")
    ggsave(p1,file=paste0("../images/epigenetic_TFBSbycommunity_OR.pdf"),height=10,width=16)

core_communities_cpg_tophits$feature<-"TFBS"
###########
load("../results/ext_communities_cpg_tophits.rdata")
#spl<-strsplit(as.character(ext_communities_cpg_tophits$filename),split=".",fixed=T)
#spl<-do.call("rbind",spl)

ext_communities_cpg_tophits<-ext_communities_cpg_tophits[which(ext_communities_cpg_tophits$collection%in%c("roadmap_epigenomics")),]
w<-which(ext_communities_cpg_tophits$logOddsRatio=="Inf")
ext_communities_cpg_tophits<-ext_communities_cpg_tophits[-w,]
ext_communities_cpg_tophits$feature<-""
g<-grep("hotspot.fdr0.01.peaks.bed",ext_communities_cpg_tophits$filename)
ext_communities_cpg_tophits$feature[g]<-"DNase hotspots"
g<-grep("DNase.macs2.narrowPeak",ext_communities_cpg_tophits$filename)
ext_communities_cpg_tophits$feature[g]<-"DNase peaks"
g<-grep("DNase.hotspot.all.peaks",ext_communities_cpg_tophits$filename)
ext_communities_cpg_tophits$feature[g]<-"DNase hotspots (all)"
w<-which(!is.na(ext_communities_cpg_tophits$antibody))
ext_communities_cpg_tophits[w,"feature"]<-"histone marks"

ext_communities_cpg_tophits<-ext_communities_cpg_tophits[which(ext_communities_cpg_tophits$userSet%in%clusters),]
ext_communities_cpg_tophits$Pvalue<-10^(-ext_communities_cpg_tophits$pValueLog)

##
w<-which(ext_communities_cpg_tophits$feature%in%"histone marks")
    m<-max(ext_communities_cpg_tophits$logOddsRatio[w])+1
    p1<-ggplot(ext_communities_cpg_tophits[w,],aes(x=antibody,y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=userSet))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
    #guides(fill = guide_legend(ncol=20))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Histone mark", fill="community")
      
    ggsave(p1,file=paste0("../images/epigenetic_histonemarks_OR.pdf"),height=10,width=16)

  p1<-ggplot(ext_communities_cpg_tophits[w,],aes(x=factor(userSet),y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=antibody))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
      guides(fill = guide_legend(nrow=2))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Community",fill="Histone mark")
    ggsave(p1,file=paste0("../images/epigenetic_histonemarksbycommunity_OR.pdf"),height=10,width=16)

##
w<-which(ext_communities_cpg_tophits$feature%in%"DNase peaks")
    m<-max(ext_communities_cpg_tophits$logOddsRatio[w])+1
    p1<-ggplot(ext_communities_cpg_tophits[w,],aes(x=antibody,y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=userSet))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
    #guides(fill = guide_legend(ncol=20))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="DNAse peak", fill="Community")
      
    ggsave(p1,file=paste0("../images/epigenetic_dnasepeaks_OR.pdf"),height=10,width=16)

  p1<-ggplot(ext_communities_cpg_tophits[w,],aes(x=factor(userSet),y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=feature))+
      #facet_wrap(~cis_snp,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
      guides(fill = guide_legend(nrow=2))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="DNAse peak",fill="Community")
    ggsave(p1,file=paste0("../images/epigenetic_dnasepeaksbycommunity_OR.pdf"),height=10,width=16)

#####
res<-rbind(ext_communities_cpg_tophits,core_communities_cpg_tophits)
res<-res[res$feature%in%c("DNase peaks","TFBS","histone marks"),]

res$userSet<-factor(res$userSet,levels=clusters)
    
    #m<-max(res$logOddsRatio)+1
    p1<-ggplot(res,aes(x=antibody,y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=userSet))+
      facet_wrap(~feature,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
    #guides(fill = guide_legend(ncol=20))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Feature", fill="community")
      
    ggsave(p1,file=paste0("../images/epigenetic_features_OR.pdf"),height=10,width=16)

  p1<-ggplot(res,aes(x=factor(userSet),y=logOddsRatio,size=-log10(Pvalue)))+
      geom_hline(yintercept=1,col="black",linetype="dotted")+
      geom_point(aes(color=antibody))+
      facet_wrap(~feature,scale="free_y",ncol=1)+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
      scale_size(range=c(1,4))+
      #scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE")) +
      guides(fill = guide_legend(nrow=2))+
      #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,10,20),limits=c(0.2,round(m))) +
      
      scale_fill_brewer(type="qual") +
      theme(legend.text=element_text(size=12)) +
    labs(y="Odds ratio (log scale)",x="Histone mark",fill="TFBS")
    ggsave(p1,file=paste0("../images/epigenetic_featuresbycommunity_OR.pdf"),height=10,width=16)


