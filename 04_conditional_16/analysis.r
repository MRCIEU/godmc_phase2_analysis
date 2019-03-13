library(tidyverse)
library(ggthemes)
library(meffil)
library(dplyr)
library(gridExtra)
library(stringr)
library(data.table)

# Conditional 

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_conditional.rdata") #1321186
#sig<-conditional[which(conditional$cis==TRUE & conditional$p<1e-8 | conditional$cis==FALSE & conditional$p<1e-14),]#1321186
sig<-conditional[which(conditional$cis==TRUE & conditional$pJ<1e-8 | conditional$cis==FALSE & conditional$pJ<1e-14),] #1058292
sig_cis<-conditional[which(conditional$cis==TRUE & conditional$pJ<1e-8),]
table(sig$cis)
nrow(sig)
nrow(sig_cis)

## Number of independent SNPs per cis and trans

clump_counts <- dplyr::group_by(sig, cpg, cis) %>%
	dplyr::summarise(numbersnps=n(),samplesize=mean(n), max_samplesize=max(n),min_samplesize=min(n))

length(which(clump_counts$samplesize>27750))
#5112

cis_counts<-clump_counts[clump_counts$cis==TRUE,]
trans_counts<-clump_counts[clump_counts$cis!=TRUE,]
median(cis_counts$numbersnps) #3
median(trans_counts$numbersnps) #1
median(clump_counts$numbersnps) #2
min(cis_counts$numbersnps) #1
min(trans_counts$numbersnps) #1
min(clump_counts$numbersnps) #1
max(cis_counts$numbersnps) #98
max(trans_counts$numbersnps) #51
max(clump_counts$numbersnps) #98

clump_counts[which(clump_counts$numbersnps>75),]
# A tibble: 2 x 4
# Groups:   cpg [2]
#  cpg        cis   numbersnps data 
#  <chr>      <lgl>      <int> <chr>
#1 cg13601595 TRUE          98 godmc
#2 cg16664915 TRUE          82 godmc


p1 <- ggplot(clump_counts, aes(x=numbersnps)) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from conditional analysis (p < 1e-8; 1e-14)", y="mQTL per DNA methylation site")
ggsave(p1, file="./images/conditional_counts_cpg.pdf", width=7, height=7)

p1 <- ggplot(clump_counts, aes(x=samplesize,y=numbersnps)) +
geom_point() +
stat_smooth(method="lm",col="red") +
labs(x="average samplesize",y="number of independent cis SNPs per DNA methylation site") +
ylim(0,100)
ggsave(plot=p1, file="./images/clump_counts_cpg_samplesize.pdf", width=7, height=7)

df<-data.frame(clump_counts)
w1<-which(df$numbersnps<5)
w2<-which(df$numbersnps>=5)
df$cat<-df$numbersnps
df$cat[w1]<-"<5 snps"
df$cat[w2]<-">=5 snps"

summary(lm(clump_counts$samplesize~clump_counts$numbersnps))
#Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             24002.473     10.281  2334.6   <2e-16 ***
#clump_counts$numbersnps  -162.191      1.218  -133.1   <2e-16 ***

p1<-ggplot(df, aes(samplesize, colour = cat)) +
geom_density()
ggsave(p1, file="./images/conditional_counts_cpg_samplesize_density.pdf", width=7, height=7)

p1<-ggplot(df, aes(max_samplesize, colour = cat)) +
geom_density()
ggsave(p1, file="./images/conditional_counts_cpg_maxsamplesize_density.pdf", width=7, height=7)

clump_counts_cis<-clump_counts[clump_counts$cis=="TRUE",]
df<-data.frame(table(clump_counts_cis$numbersnps),data="GoDMC 1000G (n=27,750)")
w<-which(as.numeric(as.character(df$Var1))<26)
s<-sum(df[-w,"Freq"])
df<-rbind(df[w,],data.frame(Var1=">25",Freq=s,data="GoDMC 1000G (n=27,750)"))
sum(df$Freq)
#[1] 181467
df$perc<-df$Freq/sum(df$Freq)

clump_counts_trans<-clump_counts[clump_counts$cis=="FALSE",]
df_trans<-data.frame(table(clump_counts_trans$numbersnps),data="GoDMC 1000G (n=27,750)")
w<-which(as.numeric(as.character(df_trans$Var1))<26)
s<-sum(df_trans[-w,"Freq"])
df_trans<-rbind(df_trans[w,],data.frame(Var1=">25",Freq=s,data="GoDMC 1000G (n=27,750)"))
sum(df_trans$Freq)
#[1] 181467
df_trans$perc<-df_trans$Freq/sum(df_trans$Freq)


#ARIES
y<-meffil.get.features("450k")
f7<-read.table("~/repo/godmc_phase2_analysis/03_clumping_16/F7.ALL.M.tab",he=T)
m<-match(f7$gene,y$name)
f7<-data.frame(cpgchr=y[m,c("chromosome")],cpgpos=y[m,c("position")],f7)
f7$cpgchr<-gsub("chr","",f7$cpgchr)
f7$cpgchr<-gsub("X","23",f7$cpgchr)
f7$cpgchr<-gsub("Y","24",f7$cpgchr)

bim<-read.table("~/repo/godmc_phase2_analysis/03_clumping_16/ariesmqtlsnps.bim")
cis_radius <- 1000000

m<-match(f7$SNP,bim[,2])
f7<-data.frame(snpchr=bim[m,1],snppos=bim[m,4],f7)

f7$cis <- FALSE
f7$cis[f7$snpchr == f7$cpgchr & (abs(f7$snppos - f7$cpgpos) <= cis_radius)] <- TRUE
sigf7 <- subset(f7, (cis & p.value < 1e-8))
clump_countsf7 <- dplyr::group_by(sigf7, gene, cis) %>%
	dplyr::summarise(numbersnps=n(),data="aries_f7")

df7<-data.frame(table(clump_countsf7$numbersnps),data="ARIES childhood 1000G (n=834)")
df7$perc<-df7$Freq/sum(df7$Freq)

#bonder
bonder<-read.table("~/repo/godmc_phase2_analysis/03_clumping_16/Bonder.txt",sep="\t",he=T)
bonder<-data.frame(Var1=bonder$noSNPs,Freq=bonder$Bonder,data="Bonder et al. GoNL (n=3,841)")
bonder$perc<-bonder$Freq/sum(bonder$Freq)

##
load("~/repo/godmc_phase2_analysis/results/16/16_conditional_3cohorts.rdata") #3984 samples
max(conditional$n) #7666.97
sig_3<-conditional[which(conditional$pJ<1e-8),] #10769

clump_counts_sub <- dplyr::group_by(sig_3, cpg) %>%
	dplyr::summarise(numbersnps=n(),data="godmc_subset")

df_sub<-data.frame(table(clump_counts_sub$numbersnps),data="GoDMC 1000G (n=3,984)")
w<-which(as.numeric(as.character(df_sub$Var1))<26)
s<-sum(df_sub[-w,"Freq"])
df_sub<-rbind(df_sub[w,],data.frame(Var1=">25",Freq=s,data="GoDMC 1000G (n=3,984)"))
sum(df_sub$Freq)
#[1] 181467
df_sub$perc<-df_sub$Freq/sum(df_sub$Freq)

##
load("~/repo/godmc_phase2_analysis/results/16/16_conditional_hrc_cis.rdata") #only cis
max(conditional.cis$n) #47140.9
sig_hrc<-conditional.cis[which(conditional.cis$pJ<1e-8),] #723524
clump_counts_hrc <- dplyr::group_by(sig_hrc, cpg) %>%
	dplyr::summarise(numbersnps=n(),data="godmc_hrc")


df_hrc<-data.frame(table(clump_counts_hrc$numbersnps),data="GoDMC HRC (n=27,750)")
w<-which(as.numeric(as.character(df_hrc$Var1))<26)
s<-sum(df_hrc[-w,"Freq"])
df_hrc<-rbind(df_hrc[w,],data.frame(Var1=">25",Freq=s,data="GoDMC HRC (n=27,750)"))
sum(df_hrc$Freq)
#[1] 181467
df_hrc$perc<-df_hrc$Freq/sum(df_hrc$Freq)
##

df_comb<-rbind(df,df7,bonder,df_sub,df_hrc)
df_comb$data<-factor(df_comb$data, levels = c("GoDMC 1000G (n=27,750)","GoDMC HRC (n=27,750)","GoDMC 1000G (n=3,984)","Bonder et al. GoNL (n=3,841)","ARIES childhood 1000G (n=834)"))

pdf("./images/conditional_counts_cpg_cis_godmc_vs_aries_bonder_subset_barplot.pdf", width=10, height=7)
p1<-ggplot(df_comb, aes(x=as.factor(Var1),y=as.numeric(perc),fill=data)) +
geom_bar(stat="identity",position="dodge") +
scale_fill_brewer(type="qual")+
labs(x="Number of SNPs per DNA methylation site",y="Fraction of DNA methylation sites with a cis association",fill="study")
#ggsave(p1, file="./images/conditional_counts_cpg_cis_godmc_vs_aries_bonder_barplot.pdf", device="pdf",width=7, height=7)
print(p1)
dev.off()


## Number of hits per SNP


sig <- subset(conditional, (cis & p < 1e-8) | (!cis & p < 1e-14))
clump_counts <- dplyr::group_by(sig, cpg, cis) %>%
	dplyr::summarise(numbersnps=n(),samplesize=mean(n))

p1 <- ggplot(clump_counts, aes(x=numbersnps)) +
geom_bar(position="dodge", aes(fill=cis)) +
labs(x="Independent hits from conditional analysis (p < 1e-8; 1e-14)", y="mQTLs per CpG")
ggsave(p1, file="./images/conditional_counts_cpg.pdf", width=7, height=7)

clump_counts_snp <- dplyr::group_by(sig, snp, cis) %>%
	dplyr::summarise(n=n()) %>%
	dplyr::group_by(n, cis) %>%
	dplyr::summarise(count=n())

p2 <- ggplot(filter(clump_counts_snp, n < 100 & n > 5), aes(x=as.factor(n), y=count)) +
geom_bar(position="dodge", aes(fill=cis), stat="identity") +
labs(x="Independent hits from conditional analysis (p < 1e-8; 1e-14)", y="mQTLs per SNP")
ggsave(p2, file="./images/conditional_counts_snp.pdf", width=7, height=7)



cl_counts <- dplyr::group_by(clumped2, cpg, cis) %>%
	dplyr::summarise(n=n())

#In the conditional analysis there is one probe with 113 cis associations (cg13601595) and 160 probes with more than 50 cisassociations. 
#In clumped there are no probes with more than 50 associations and only 29 probes with more than 4 associations.
