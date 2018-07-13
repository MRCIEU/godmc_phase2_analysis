library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/enrichments/mqtl_pchic.rdata")
#snp in bait; cpg in oe
snp_in_promoter<-unique(data$snp_in_promoter)

#cpg in bait; snp in oe
cpg_in_promoter<-unique(data$cpg_in_promoter)

ann[ann$Name%in%c(snp_in_promoter$CpG[1:3]),1:5]


load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/05_cis-trans-networks/results/communities.rdata")
c8<-communities[communities$cluster==8,]
c8[which(c8$cpg%in%cpg_in_promoter$CpG),]
c8[which(c8$snp%in%snp_in_promoter$SNP),]
c8$code<-paste(c8$cpg,c8$snp)

hic<-unique(c(snp_in_promoter$code,cpg_in_promoter$code))
length(which(c8$code%in%hic))

nrow(communities) #13651
spl<-do.call("rbind",strsplit(communities$snp,split=":"))
w<-which(spl[,1]=="chr6"&spl[,2]>29570005&spl[,2]<33377657)
communities<-communities[-w,]
nrow(communities) #12405

communities<-data.frame(communities)
#clusters<-c("8","31","12")
df<-data.frame(table(communities$cluster)) #2316
clusters<-df[which(df$Freq>10),1] #165
communities$code<-paste(communities$cpg,communities$snp)

hic.res<-data.frame()
for (i in 1:length(clusters)){
cat(i,"\n")

cl<-unique(communities[which(communities$cluster==clusters[i]),6])
hic.res[i,1]<-clusters[i]
hic.res[i,2]<-length(which(cl%in%hic))/length(unique(cl))
}

l<-length(which(unique(communities$code)%in%hic))/length(unique(communities$code))
#0.003

hic.res[which(hic.res[,2]>l),]

hic.res[which(hic.res[,2]>l),]
#     V1          V2
#6     8 0.013888889
#11   14 0.023255814
#26   32 0.071428571
#29   36 0.007633588
#52   68 0.052631579
#70  105 0.090909091
#86  134 0.030303030
#98  162 0.125000000
#103 182 0.013157895
#118 214 0.043478261
#119 215 0.068181818

hic.res[which(hic.res[,2]>0.05),1]
load("../results/KEGGenrichments_fdr0.05.rdata")
df.out[df.out$cluster%in%hic,]
load("../results/goenrichments_fdr0.05.rdata")

load("../results/core_communities_cpg_tophits.rdata")
data.frame(core_communities_cpg_tophits[which(core_communities_cpg_tophits$userSet%in%hic),])

load("../results/ext_communities_cpg_tophits.rdata")
tiss<-read.table("../data/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
spl<-do.call("rbind",strsplit(ext_communities_cpg_tophits$filename,split="-"))
ext_communities_cpg_tophits$dataSource<-spl[,1]
m<-match(ext_communities_cpg_tophits$dataSource,as.character(tiss$ID),)
ext_communities_cpg_tophits$tissue<-tiss[m,"Tissue"]

data.frame(ext_communities_cpg_tophits[which(ext_communities_cpg_tophits$userSet%in%hic),])
#  userSet dbSet          collection pValueLog logOddsRatio support rnkPV rnkLO
#1      32   466 roadmap_epigenomics  4.975731    11.629650      13     1     4
#2      32  1149 roadmap_epigenomics  4.351454     8.496724      11     2     9
#3      32   480 roadmap_epigenomics  4.163346     9.294153      13     3     7
#4      32   608 roadmap_epigenomics  4.005042    10.874158      14     4     6
#  rnkSup maxRnk meanRnk    b c    d                 description cellType
#1      7      7    4.00 1382 3 3711 roadmap_epigenomics H3K4me1     <NA>
#2     14     14    8.33 1047 5 4046 roadmap_epigenomics H3K4me1     <NA>
#3      7      7    5.67 1619 3 3474 roadmap_epigenomics H3K4me1     <NA>
#4      5      6    5.00 1994 2 3099 roadmap_epigenomics H3K4me1     <NA>
#    tissue antibody treatment dataSource                filename   size
#1    BLOOD  H3K4me1      <NA>       E030 E030-H3K4me1.narrowPeak 116781
#2 VASCULAR  H3K4me1      <NA>       E122 E122-H3K4me1.narrowPeak 142297
#3    BLOOD  H3K4me1      <NA>       E032 E032-H3K4me1.narrowPeak 167667
#4    BLOOD  H3K4me1      <NA>       E051 E051-H3K4me1.narrowPeak 223023
#        fdr       fdr2
#1 0.0805206 0.01311264
#2 0.2068536 0.02760182
#3 0.2805772 0.02837623
#4 0.3629311 0.03064220

#snp_in_promoter <- list(table=dat, test=fisher, unique_snps=length(unique(snp_in_promoter$SNP)), percentage=length(unique(snp_in_promoter$SNP)) / length(unique(bait_pchic_snp$SNP)[!unique(bait_pchic_snp$SNP) %in% unique(snp_in_promoter$SNP)]))
#cpg_in_promoter <- list(table=dat, test=fisher, unique_cpgs=length(unique(cpg_in_promoter$CpG)), percentage=length(unique(cpg_in_promoter$CpG)) / length(unique(bait_pchic_cpg$CpG)[!unique(bait_pchic_cpg$CpG) %in% unique(cpg_in_promoter$CpG)]))

