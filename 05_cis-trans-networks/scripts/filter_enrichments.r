library(ggplot2)
library(dplyr)

load("../results/communities.rdata")
ncomm<-communities %>%
  group_by(cluster) %>%
  summarize(n=n())

length(which(ncomm$n>2))
#1331
#load("../results/gwas_clusters.rdata")
load("../results/gwas_clusters_nochr6.rdata")

info <- read.csv("../../data/gwas/00info.csv")
info <- data.frame(fn=gsub(".txt.gz", "", info$newfile), id=info$id)

labels <- subset(res, !duplicated(i), select=c(i, fn))
labels$fn <- gsub("../../data/gwas/", "", labels$fn)
labels$fn <- gsub(".txt.gz", "", labels$fn)
dat <- merge(dat, labels, by.x="id", by.y="i")
dat$fisher <- unlist(dat$fisher)
dat$fdr <- p.adjust(dat$fisher, "fdr")
dat$bonferroni <- p.adjust(dat$fisher, "bonferroni")
dat$bonferroni3 <- p.adjust(dat$binom4, "bonferroni")
dat$label <- gsub("disease__", "", dat$fn)
dat$label <- gsub("risk_factor__", "", dat$label)
dat$label <- gsub("_", " ", dat$label)


dat <- merge(dat, info, by="fn")
load("../data/outcomes.RData")
ao <- subset(ao, select=c(id, subcategory))
dat <- merge(dat, ao, by.x="id.y", by.y="id",all.x=T)
#dat <- merge(dat, ao, by.x="id.y", by.y="id")
unique(dat[which(is.na(dat$subcategory)),c("id.y","fn","label")])

dat$subcategory[dat$subcategory=="Hemodynamic"] <- "Haemotological"
dat$subcategory[dat$subcategory=="Immune system"] <- "Autoimmune / inflammatory"
dat$subcategory[dat$subcategory=="Diabetes"] <- "Glycemic"
dat$subcategory[dat$subcategory=="Biomarker"] <- "Other"
dat$subcategory[dat$subcategory=="Protein"] <- "Other"
dat$subcategory[dat$subcategory=="Reproductive aging"] <- "Aging"
dat$subcategory[dat$subcategory=="Lung disease"] <- "Other"
dat$subcategory[dat$subcategory=="Autoimmune / inflammatory"] <- "Immune"
dat$subcategory[dat$subcategory=="Psychiatric / neurological"] <- "Neurological"

gwa<-group_by(dat, clust) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% as.data.frame

g<-grep("metabolites__",dat$fn)

dat2<-dat[-g,]
dat2<-dat2[which(dat2$nsnp>5),]

temp <- group_by(dat, id.y) %>%
summarise(
	nsnp=mean(nsnp),
	nsig=sum(fisher < 0.05),
	fdr=sum(fdr < 0.05),
	bonferroni=sum(bonferroni3 < 0.05),
	minp=min(min_p, na.rm=TRUE)
) %>% filter(bonferroni > 0)

g<-grepl("metabolites__", dat$fn) 
dat2<-dat[g==FALSE,] #1843
dat2<-dat2[which(dat2$nsnp>=5),] #1843
dat2<-dat2[which(dat2$id.y %in% temp$id.y),] #302

min(dat2$binom4)
#[1] 3.048592e-39
o<-order(dat2$binom4)
dat2<-dat2[o,]

g<-grepl("metabolites__", dat$fn) 
dat3<-dat[g==FALSE,] #1843
dat3<-dat3[which(dat3$nsnp>3),] #1843

sign<-0.05/nrow(dat)
#[1] 5.156766e-06

dat4<-dat3[dat3$binom4<sign,] #43
#table(dat4$clust)

# 3  8 12 14 15 21 23 26 31 35 36 43 51 
# 1 10  1  2  1 14  1  1  1  8  1  1  1 

#cluster,hla,all snps 
#3,6,6
#8,32,33
#12,0,6
#14,23,23
#15,6,6
#21,17,22
#23,2,10
#26,7,8
#31,0,10
#35,15,15
#36,4,13
#43,8,10
#51,5,6

n<-names(table(dat4$clust))

for (i in 1:length(n)){
cat(nrow(unique(communities[communities$cluster==n[i],"snp"])),"\n")
}

clust<-unique(dat4$clust)

p1<-ggplot(dat3, aes(x=label, y=-log10(binom4))) +
geom_point(aes(colour=clust, size=nsnp)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
geom_hline(yintercept=-log10(0.05/nrow(dat)), linetype="dotted") +
labs(x="", y="Enrichment", size="Number\nof SNPs in\ncommunity") +
scale_colour_continuous(guide=FALSE) +
facet_grid(. ~ subcategory, scale="free", space="free") +
theme(legend.position="bottom", strip.text=element_text(angle=90, size=10), axis.text.x=element_text(size = 10))
ggsave(p1,file="../images/gwas_clusters_full_test.pdf", width=12, height=10)


min(dat3$binom4)
#[1] 3.048592e-39
o<-order(as.numeric(dat3$binom4))
dat3<-dat3[o,]

#id.y = ao$id

load("../results/core_communities_cpg_tophits.rdata")
tiss<-read.table("../data/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
spl<-do.call("rbind",strsplit(core_communities_cpg_tophits$filename,split="-"))
core_communities_cpg_tophits$dataSource<-spl[,1]
m<-match(core_communities_cpg_tophits$dataSource,as.character(tiss$ID),)
core_communities_cpg_tophits$tissue<-tiss[m,"Tissue"]

c8<-which(core_communities_cpg_tophits$userSet%in%clust)
c8<-core_communities_cpg_tophits[c8,]
table(c8$userSet,c8$antibody)
    
#     ATF3 BDP1 BRF1 EZH2_(39875) RPC155
#  8     1    1    1            0      1
#  12    0    0    0            1      0

c8_2<-c8[which(c8$fdr<0.05),]
table(c8_2$userSet,c8_2$antibody)
#    BDP1 RPC155
#  8    1      1

#BDP1:Subunit Of RNA Polymerase III Transcription Initiation Factor IIIB
#RPC155:RNA Polymerase III Subunit A

# A tibble: 5 x 24
# Groups:   userSet [2]
#  userSet dbSet  collection pValueLog logOddsRatio support rnkPV rnkLO rnkSup
#    <int> <int>       <chr>     <dbl>        <dbl>   <int> <int> <int>  <int>
#1       8   985 encode_tfbs  8.575652    146.91644       5     1     7     77
#2       8   925 encode_tfbs  5.861692    262.80446       3     2     4    184
#3      12   461 encode_tfbs  4.302928     12.01578       6     2    17     24
#4       8   927 encode_tfbs  3.976064    249.03860       2     3     5    293
#5       8   500 encode_tfbs  3.568651     15.44613       4     6    33    108
# ... with 15 more variables: maxRnk <int>, meanRnk <dbl>, b <int>, c <int>,
#   d <int>, description <chr>, cellType <chr>, tissue <chr>, antibody <chr>,
#   treatment <chr>, dataSource <chr>, filename <chr>, size <dbl>, fdr <dbl>,
#   fdr2 <dbl>
#> data.frame(c8)
#  userSet dbSet  collection pValueLog logOddsRatio support rnkPV rnkLO rnkSup
#1       8   985 encode_tfbs  8.575652    146.91644       5     1     7     77
#2       8   925 encode_tfbs  5.861692    262.80446       3     2     4    184
#3      12   461 encode_tfbs  4.302928     12.01578       6     2    17     24
#4       8   927 encode_tfbs  3.976064    249.03860       2     3     5    293
#5       8   500 encode_tfbs  3.568651     15.44613       4     6    33    108
#  maxRnk meanRnk   b  c    d               description cellType tissue
#1     77    28.3  10 17 5077          ChIP K562 RPC155     K562   <NA>
#2    184    63.3   3 19 5084            ChIP K562 BDP1     K562   <NA>
#3     24    14.3 188 13 4902 ChIP H1-hESC EZH2_(39875)  H1-hESC   <NA>
#4    293   100.0   2 20 5085            ChIP K562 BRF1     K562   <NA>
#5    108    49.0  72 18 5015            ChIP A549 ATF3     A549   <NA>
#      antibody    treatment dataSource
#1       RPC155         None       <NA>
#2         BDP1         None       <NA>
#3 EZH2_(39875)         None       <NA>
#4         BRF1         None       <NA>
#5         ATF3 EtOH_0.02pct       <NA>
#                                                   filename size          fdr
#1             wgEncodeAwgTfbsSydhK562Rpc155UniPk.narrowPeak 1310 3.455065e-05
#2               wgEncodeAwgTfbsSydhK562Bdp1UniPk.narrowPeak  570 1.125022e-02
#3       wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak 6370 2.870163e-01
#4               wgEncodeAwgTfbsSydhK562Brf1UniPk.narrowPeak  221 4.982138e-01
#5 wgEncodeAwgTfbsHaibA549Atf3V0422111Etoh02UniPk.narrowPeak 6580 9.131560e-01
#          fdr2
#1 1.830491e-06
#2 4.736935e-04
#3 3.429978e-02
#4 2.426800e-02
#5 4.650593e-02


load("../results/ext_communities_cpg_tophits.rdata")
tiss<-read.table("../data/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
spl<-do.call("rbind",strsplit(ext_communities_cpg_tophits$filename,split="-"))
ext_communities_cpg_tophits$dataSource<-spl[,1]
m<-match(ext_communities_cpg_tophits$dataSource,as.character(tiss$ID),)
ext_communities_cpg_tophits$tissue<-tiss[m,"Tissue"]


e8<-which(ext_communities_cpg_tophits$userSet%in%clust)
e8<-ext_communities_cpg_tophits[e8,]
table(e8$userSet,e8$antibody)
    
#     H2AK5ac H2A.Z H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K56ac
#  3        0     0       0        0       60       5       0       0       0
#  8        0     0       0        0        0       0       0       0       0
#  12       1     0       2       65        0      35       0       1       1
#  21       0     1      10        0        0      42       4       7       0
#  23       0     0       0        0        0       0       0       0       0
#  31       0     0       0        0       17       0       0       0       0
    
#     H3K79me2 H3K9ac H3K9me3 H4K20me1
#  3         1      0       0        0
#  8         0      0       0        0
#  12        5      0       1        3
#  21        0      9       0        0
#  23        0      0       0        0
#  31        0      0       0        0

e8_2<-e8[which(e8$fdr<0.05),]
table(e8_2$userSet,e8_2$antibody)
    
#     H3K27me3 H3K36me3 H3K4me1
#  3         0       13       0
#  12        6        0       0
#  21        0        0       4
#  31        0        1       0

#H3K36me3:linked to active transcription
dat4[dat4$clust==3,c("label","binom4")]
#                    label       binom4
#5975 rheumatoid arthritis 1.995504e-08

dat4[dat4$clust==21,c("label","binom4")]
#                            label       binom4
#59                  schizophrenia 1.254968e-10
#264                        height 8.378419e-17
#380     haemoglobin concentration 8.362276e-08
#446                     psoriasis 4.996000e-12
#473  systemic lupus erythematosus 1.254968e-10
#505               hdl cholesterol 1.995504e-08
#529             total cholesterol 5.995000e-15
#541                 triglycerides 5.995000e-15
#5928     pgc crossdisorder traits 5.586012e-14
#5956         rheumatoid arthritis 2.413552e-23
#8640                  lung cancer 8.362276e-08
#8648    squamous cell lung cancer 8.362276e-08
#8896             age at menopause 6.977628e-11
#9587         mothers age at death 8.068712e-07


load("../results/s25_communities_cpg_tophits.rdata")
spl<-do.call("rbind",strsplit(s25_communities_cpg_tophits$filename,split="_"))
s25_communities_cpg_tophits$dataSource<-spl[,1]
s25_communities_cpg_tophits$seg_code<-gsub(".bed","",spl[,5])

state<-read.table("../data/roadmap_states.txt",he=T)
state$STATE<-paste0("E",state$STATE)
m<-match(s25_communities_cpg_tophits$seg_code,state$STATE)
s25_communities_cpg_tophits$seg_explanation<-state[m,"NO."]

tiss<-read.table("../data/jul2013.roadmapData_tissues.txt",sep="\t",he=T)
m<-match(s25_communities_cpg_tophits$dataSource,as.character(tiss$ID),)
s25_communities_cpg_tophits$tissue<-tiss[m,"Tissue"]


s8<-which(s25_communities_cpg_tophits$userSet%in%clust)
s8<-s25_communities_cpg_tophits[s8,]
#PromU=Promoter Upstream TSS
table(data.frame(s8[s8$tissue=="BLOOD","userSet"]))
# 3 12 31 
#11  6  9 
table(s8$userSet,s8$seg_explanation)
    
#     DNase EnhA1 EnhA2 EnhAc EnhAF EnhW1 EnhW2 Het PromBiv PromD1 PromD2 PromP
#  3      0     0     0     3     0     0     0   0       0      0      0     0
#  8      0     0     0     0     0     0     0   0       0      0      0     0
#  12     0     0     0     0     0     0     0   0       4      1     13     0
#  21     0     1     0     0     0     0     0   0       0      0      0     0
#  23     0     0     0     0     0     0     0   0       0      0      0     0
#  31     0     0     0     0     0     0     0   0       0      0      0     0
#  36     0     0     0     0     0     0     1   0       0      0      0     0
    
 #    PromU Quies ReprPC TssA Tx Tx3prime Tx5prime TxEnh3prime TxEnh5prime
 # 3      0     0      0    0 36       29        0           4          14
 # 8      5     0      0    0  0        0        0           0           0
 # 12     0     0      0    0  0        0        0           0           0
 # 21     0     0      0    0  0        0        0           0           0
 # 23     0     0      0    2  0        0        0           0           0
 # 31     0     0      0    0  0       35        0           0           0
 # 36     0     0      0    0  0        0        0           0           0
    
 #    TxEnhW TxReg TxWk ZNF/Rpts
 # 3       3     0    0        0
 # 8       0     0    0        0
 # 12      0    19    0        0
 # 21      0     0    0        0
 # 23      0     0    0        0
 # 31      0     0   17        0
 # 36      0     0    0        0

s8_2<-s8[which(s8$fdr<0.05),]
table(s8_2$userSet,s8_2$seg_explanation)
    
#     DNase EnhA1 EnhA2 EnhAc EnhAF EnhW1 EnhW2 Het PromBiv PromD1 PromD2 PromP
#  3      0     0     0     0     0     0     0   0       0      0      0     0
#  8      0     0     0     0     0     0     0   0       0      0      0     0
#  12     0     0     0     0     0     0     0   0       0      0      0     0
#  21     0     1     0     0     0     0     0   0       0      0      0     0
#  31     0     0     0     0     0     0     0   0       0      0      0     0
    
#     PromU Quies ReprPC TssA Tx Tx3prime Tx5prime TxEnh3prime TxEnh5prime
#  3      0     0      0    0  2        0        0           0           0
#  8      2     0      0    0  0        0        0           0           0
#  12     0     0      0    0  0        0        0           0           0
#  21     0     0      0    0  0        0        0           0           0
#  31     0     0      0    0  0        1        0           0           0
    
#     TxEnhW TxReg TxWk ZNF/Rpts
#  3       0     0    0        0
#  8       0     0    0        0
#  12      0     2    0        0
#  21      0     0    0        0
#  31      0     0    0        0


data.frame(table(s25_communities_cpg_tophits$seg_explanation))
#          Var1 Freq
#1        DNase  149
#2        EnhA1   71
#3        EnhA2   66
#4        EnhAc   14
#5        EnhAF   22
#6        EnhW1   47
#7        EnhW2   53
#8          Het  248
#9      PromBiv   63
#10      PromD1    8
#11      PromD2   30
#12       PromP    2
#13       PromU    5
#14       Quies  187
#15      ReprPC  361
#16        TssA  197
#17          Tx   42
#18    Tx3prime   64
#19    Tx5prime    1
#20 TxEnh3prime   29
#21 TxEnh5prime   83
#22      TxEnhW   40
#23       TxReg  201
#24        TxWk   18
#25    ZNF/Rpts    0

load("../results/communities.rdata")
nrow(communities) #13651
spl<-do.call("rbind",strsplit(communities$snp,split=":"))
w<-which(spl[,1]=="chr6"&spl[,2]>29570005&spl[,2]<33377657)
communities<-communities[-w,]
nrow(communities) #12405

communities<-data.frame(communities)
clusters<-c("8","31","12")
for (i in 1:length(clusters)){
sigCpGs<-unique(as.character(communities[which(communities$cluster==clusters[i]),1]))
g8 <- gometh(sig.cpg=sigCpGs, all.cpg=unique(communities$cpg))
print(g8[g8$FDR<0.05,])}

#g8[g8$FDR<0.05,]
#                                               Term Ont  N DE         P.DE
#GO:0006323                            DNA packaging  BP 49  5 2.381259e-05
#GO:0006333        chromatin assembly or disassembly  BP 47  5 1.295659e-05
#GO:0006334                      nucleosome assembly  BP 37  5 4.168682e-06
#GO:0031497                       chromatin assembly  BP 40  5 6.168938e-06
#GO:0034728                  nucleosome organization  BP 44  5 1.352560e-05
#GO:0065004             protein-DNA complex assembly  BP 50  5 1.137304e-05
#GO:0071824 protein-DNA complex subunit organization  BP 57  5 2.995467e-05
#GO:0000786                               nucleosome  CC 30  5 5.566210e-07
#GO:0032993                      protein-DNA complex  CC 49  6 1.145897e-06
#GO:0044815                    DNA packaging complex  CC 31  5 5.979691e-07
#                   FDR
#GO:0006323 0.037274634
#GO:0006333 0.023818578
#GO:0006334 0.014682099
#GO:0031497 0.017381600
#GO:0034728 0.023818578
#GO:0065004 0.023818578
#GO:0071824 0.042200138
#GO:0000786 0.004212094
#GO:0032993 0.005381133
#GO:0044815 0.004212094

#[1] Term Ont  N    DE   P.DE FDR 
#<0 rows> (or 0-length row.names)
#[1] Term Ont  N    DE   P.DE FDR 
#<0 rows> (or 0-length row.names)





