library(tidyverse)
load("../results/enrichments/snpcontrolsets.rdata")

r<-read_delim("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/iHS_CEU.whole_genome.pvalues.gz",delim=" ")
names(r)[4:5]<-c("iHS_score","iHS_pval")
#[1] 7291078       5

w1<-(which(r$iHS_score>2))
length(w1)
#[1] 378265

#378265/7291078
#[1] 0.05188053

test<-r[w1,]
mean(test$iHS_pval)
#[1] 1.683528
test<-r[-w1,]
mean(test$iHS_pval)
#[1] 0.355008

w1<-(which(r$iHS_pval>2))
test<-r[w1,]

length(w1)
#[1] 66426
min(test$iHS_score)
#[1] 2.920335
max(test$iHS_score)
#[1] 20.79328

mean(test$iHS_score)
#[1] 3.386194
test<-r[-w1,]
mean(test$iHS_score)
#[1] 0.7155552

w1<-which(r$iHS_pval>2)
w2<-which(r$iHS_pval<2)
r$ihs[w2]<-0
r$ihs[w1]<-1

table(r$ihs)
#0       1 
#7224652   66426



r2<-read_delim("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/FstGLOB_CEU_u_YRI_u_CHB.whole_genome.pvalues.gz",delim=" ")
names(r2)[4:5]<-c("Fst_score","fst_pval")

w1<-(which(r2$fst_pval>2))
#[1] 205790
test<-r2[w1,]

min(test$Fst_score)
#[1] 0.3997
max(test$Fst_score)
#[1] 0.966

mean(test$Fst_score)
#[1] 0.4885446
test<-r2[-w1,]
mean(test$Fst_score)
#[1] 0.09290992

w1<-which(r2$fst_pval>2)
w2<-which(r2$fst_pval<2)
r2$fst[w2]<-0
r2$fst[w1]<-1

r3<-read_delim("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/XPEHH_CEU_vs_CHB.whole_genome.pvalues.gz",delim=" ")
names(r3)[4:5]<-c("xpehhchb_score","xpehhchb_pval")

w1<-(which(r3$xpehhchb_pval>2))
#[1] 64489
test<-r3[w1,]

min(test$xpehhchb_score)
#[1] 2.432219
max(test$xpehhchb_score)
#[1] 7.942488

mean(test$xpehhchb_score)
#[1] 2.89147
test<-r3[-w1,]
mean(test$xpehhchb_score)
#[1] -0.08157803

w1<-which(r3$xpehhchb_pval>2)
w2<-which(r3$xpehhchb_pval<2)
r3$xpehhchb[w2]<-0
r3$xpehhchb[w1]<-1

table(r3$xpehhchb)
#      0       1 
#7623540   64489 


r4<-read_delim("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/XPEHH_CEU_vs_YRI.whole_genome.pvalues.gz",delim=" ")
names(r4)[4:5]<-c("xpehhyri_score","xpehhyri_pval")

w1<-(which(r4$xpehhyri_pval>2))
#[1] 56181
test<-r4[w1,]

min(test$xpehhyri_score)
#[1] 2.918218
max(test$xpehhyri_score)
#[1] 7.333925

mean(test$xpehhyri_score)
#[1] 3.463855
test<-r4[-w1,]
mean(test$xpehhyri_score)
#[1] -0.1396265

w1<-which(r4$xpehhyri_pval>2)
w2<-which(r4$xpehhyri_pval<2)
r4$xpehhyri[w2]<-0
r4$xpehhyri[w1]<-1

table(r4$xpehhyri)

# 0       1 
#7592627   56181

r5<-read_delim("/panfs/panasas01/shared-godmc/1kg_reference_ph3/selection_results/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz",delim=" ")
names(r5)[4:5]<-c("sds_score","sds_pval")
w1<-(which(r5$sds_pval>2))
#[1] 44361
test<-r5[w1,]

min(test$sds_score)
#2.384468
max(test$sds_score)
#[1] 10.00023

mean(test$sds_score)
#[1] 2.835383
test<-r5[-w1,]
mean(test$sds_score)
#[1] -0.02781859

w1<-which(r5$sds_pval>2)
w2<-which(r5$sds_pval<2)
r5$sds[w2]<-0
r5$sds[w1]<-1

table(r5$sds)

m<-match(f.all$SNP,r$snpID)
m2<-match(f.all$SNP,r2$snpID)
m3<-match(f.all$SNP,r3$snpID)
m4<-match(f.all$SNP,r4$snpID)
m5<-match(f.all$SNP,r5$snpID)

f.all<-data.frame(f.all,r[m,c("iHS_score","iHS_pval","ihs")],r2[m2,c("Fst_score","fst_pval","fst")],r3[m3,c("xpehhchb_score","xpehhchb_pval","xpehhchb")],r4[m4,c("xpehhyri_score","xpehhyri_pval","xpehhyri")],r5[m5,c("sds_score","sds_pval","sds")])
save(f.all,file="../results/enrichments/snpcontrolsets_selection.rdata")
