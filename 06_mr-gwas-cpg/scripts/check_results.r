###
library(meffil)
library(data.table)

y<-meffil.get.features("450k")
smok<-read.table("../data/smok.csv",he=T,sep=",")
smok$Location..hg19.<-gsub(",","",smok$Location..hg19.)
smok<-smok[smok$P.value..Current.vs..Never.<1e-7,]

r<-read.csv("../results/trait-cpg-sig-codes.csv")
r<-r[r$nsnp>1&r$pval<1.4e-7,]
nrow(r) #194
length(unique(r$trait)) #14

r2<-r[r$msig==TRUE,]
nrow(r2) #85
length(unique(r2$trait)) #13

#y<-extract_instruments(outcomes=279)
#y$SNP
#[1] "rs3763350" "rs477515"  "rs6679677"

jia<-grep("Juvenile idiopathic arthritis",r2$trait)
r2<-r2[-jia,]
nrow(r2) #23

df<-data.frame(table(r2$trait))
df[df$Freq>0,]
#                         Var1 Freq
#1             Age at menarche    1
#3                Birth weight    1
#4             Bulimia nervosa    1
#5                      Copper    1
#7             r2[]sulin    1
#12                       Iron    1
#18         Multiple sclerosis    2
#25       Rheumatoid arthritis    3
#26 Serum cystatin C (eGFRcys)    2
#29     Transferrin Saturation    4
#30              Triglycerides    5
#35                       Zinc    1


m<-match(r2$outcome,y$name)
r2<-data.frame(y[m,c("name","chromosome","position","gene.symbol")],r2)
o<-order(r2$pval)
r2[o,]

hla<-which(r2$chromosome=="chr6"&r2$position>24570005&r2$position<38377657)
r3<-r2[-hla,]
hla<-r2[hla,]#9
length(unique(hla$trait))
length(unique(hla$name))

o<-order(r3$chromosome,r3$position)
r3<-r3[o,]

w<-which(r3$name%in%smok$Name)
r3[w,] #0

df.out<-data.frame()
for (i in (1:nrow(r3))){
chr<-gsub("chr","",r3[i,"chromosome"])
w<-which(smok$Chromosome==chr)
smok.chr<-smok[w,]
dist<-abs(as.numeric(smok.chr$Location..hg19.)-r3$position[i])
w<-which(dist<100000)
if(length(w)>0){
df<-data.frame(r3[i,"name"],smok.chr[w,])
df.out<-rbind(df.out,df)
}
}
w<-which(r3$name%in%unique(df.out[,1]))
r3.smok<-r3[w,] #12
r3.nonsmok<-r3[-w,]

r3<-merge(r3,y[,c("name","")],by.x=name,by.y=name)

##################################

load("/newhome/epzjlm/repo/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))
bim<-fread("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig")

#mQTL for smoking sites
clsmok<-clumped[which(clumped$cpg%in%smok$Name),]
length(unique(clsmok$snp)) #3414
length(unique(clsmok$cpg)) #2084

#age at menarche snps
gwa<-read.table("age_at_menarche.snps")
bim2<-bim[bim$V2%in%gwa$V1,] 

df.out<-data.frame()
for (i in (1:nrow(bim2))){
chr<-bim2$V1[i]
w<-which(clsmok$snpchr==paste0("chr",chr))
smok.chr<-clsmok[w,]
dist<-abs(as.numeric(smok.chr$snppos)-bim2$V4[i])
w<-which(dist<1000000)
if(length(w)>0){
df<-data.frame(bim2[i,],smok.chr[w,])
df.out<-rbind(df.out,df)
}
}

df.out<-df.out[which(df.out$Name%in%unique(clsmok$cpg)),]
nrow(bim2)#68

length(unique(df.out$V2))
#[1] 33/68 age at menarche were with 1Mb of a smoking mQTL
length(unique(df.out$V1))
# 18

gwa<-read.table("glucose_instruments.txt")
bim2<-bim[bim$V2%in%gwa$V1,]

df.out<-data.frame()
for (i in (1:nrow(bim2))){
chr<-bim2$V1[i]
w<-which(clsmok$snpchr==paste0("chr",chr))
smok.chr<-clsmok[w,]
dist<-abs(as.numeric(smok.chr$snppos)-bim2$V4[i])
w<-which(dist<1000000)
if(length(w)>0){
df<-data.frame(bim2[i,],smok.chr[w,])
df.out<-rbind(df.out,df)
}
}
nrow(bim2)
#22
length(unique(df.out$V2))
#[1] 13
length(unique(df.out$V1))
# 8

gwa<-read.table("creatinine.snps")
bim2<-bim[bim$V2%in%gwa$V1,]

df.out<-data.frame()
for (i in (1:nrow(bim2))){
chr<-bim2$V1[i]
w<-which(clsmok$snpchr==paste0("chr",chr))
smok.chr<-clsmok[w,]
dist<-abs(as.numeric(smok.chr$snppos)-bim2$V4[i])
w<-which(dist<1000000)
if(length(w)>0){
df<-data.frame(bim2[i,],smok.chr[w,])
df.out<-rbind(df.out,df)
}
}

nrow(bim2)
#46
length(unique(df.out$V2))
#31
length(unique(df.out$V1))
#16

gwa<-read.table("triglycerides.snps")
bim2<-bim[bim$V2%in%gwa$V1,]

df.out<-data.frame()
for (i in (1:nrow(bim2))){
chr<-bim2$V1[i]
w<-which(clsmok$snpchr==paste0("chr",chr))
smok.chr<-clsmok[w,]
dist<-abs(as.numeric(smok.chr$snppos)-bim2$V4[i])
w<-which(dist<1000000)
if(length(w)>0){
df<-data.frame(bim2[i,],smok.chr[w,])
df.out<-rbind(df.out,df)
}
}

nrow(bim2)
#54
length(unique(df.out$V2))
#40
length(unique(df.out$V1))
#16





