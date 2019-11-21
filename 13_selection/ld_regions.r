library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(meffil)
library(data.table)
library(ggplot2)
library(GenomicAlignments)
library(Rsamtools)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicFeatures")

#FUNCTION from MATT -really fast!
genomic.nearest <- function(positions, intervals) {
     stopifnot(all(c("chromosome","start","end") %in% colnames(intervals)))
     stopifnot(all(c("chromosome","position") %in% colnames(positions)))
 
     events <- rbind(data.frame(chromosome=intervals$chromosome,
                                position=intervals$start,
                                type="start",
                                id=1:nrow(intervals)),
                     data.frame(chromosome=intervals$chromosome,
                                position=intervals$end,
                                type="end",
                                id=1:nrow(intervals)),
                     data.frame(chromosome=positions$chromosome,
                                position=positions$position,
                                type="position",
                                id=1:nrow(positions)))
     events <- events[order(events$chromosome, events$position, decreasing=F),]
 
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
     dist.before[which(events$chromosome != events$chromosome[before])] <- NA
     dist.after <- events$position[after] - events$position
     dist.after[which(events$chromosome != events$chromosome[after])] <- NA
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

# Get clumped data


load("../results/16/16_clumped.rdata")
cl<-data.frame(clumped)
cl.out<-data.frame()
 
#
#retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
#excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
#rm<-which(retaincpg%in%excl[,1])
#14882
#retaincpg<-retaincpg[-rm]
#420509
 
#nrow(cl)
#284819

#cl<-cl[which(cl$cpg%in%retaincpg),]
#nrow(cl)
#[1] 284819

length(unique(cl$snp))
#[1] 235439

#remove snps
#flip<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/data/ref/flipped_snps.txt",he=F)
#w<-which(cl$snp%in%flip[,1])
#length(w) #295
#if(length(w)>0){
#cl<-cl[-w,]}

#indels<-read.table("/panfs/panasas01/shared-godmc/INDELs/indels_equal_seq_length.txt")
#w<-which(cl$snp%in%indels[,1]) #0
#length(w)
#if(length(w)>0){
#cl<-cl[-w,]}

#Filter on pvalue
clumped <- subset(clumped, (pval < 1e-14 & cis == FALSE) | (pval < 1e-8 & cis == TRUE ))
#nrow(clumped)
#[1] 271724
#

#data<-data.table(clumped)

#data$cis<-as.logical(data$cis)
#data$abs_Effect<-abs(data$Effect)


#df1<-data[,snp_cis:=ifelse(all(cis),"TRUE",ifelse(all(!cis),"FALSE","ambivalent")),by=c("snp")]

#df2<-data[ , (min_pval = min(pval)), by = snp]

#data<-inner_join(df1,df2)
#w<-which(names(data)%in%"V1")
#names(data)[w]<-"min_pval"
#data<-data.table(data)
#df3<-data[ , (min_Effect = min(Effect)), by = snp]
#df4<-data[ , (max_Effect = max(Effect)), by = snp]
#df<-data.table(rbind(df3,df4))

#maxAbsObs <- function(x) x[which.max(abs(x))]
#df2<-df[, lapply(.SD, maxAbsObs), by="snp"]
#data<-inner_join(data,df2)
#w<-which(names(data)%in%"V1")
#names(data)[w]<-"max_abs_Effect"

#data$mQTL<-TRUE

cl.snp<-unique(data.frame(mqtl_clumped=clumped$snp))
dim(cl.snp)
#224648

y<-meffil.get.features("450k")


#get TSS data
#distance to TSS

tss<-read.table("/panfs/panasas01/shared-godmc/gencode/gencode.v27lift37.annotation.gtf.gz",skip=5,sep="\t")
names(tss)[1:5]<-c("chromosome","source","type","start","end")

table(tss$type)


#           CDS           exon           gene Selenocysteine    start_codon 
#        713931        1206031          60461            119          83869 
#    stop_codon     transcript            UTR 
#         75703         202697         286464 

 
tss<-tss[tss$type=="transcript",]

g<-grep("gene_type protein_coding",tss$V9)
tss<-tss[g,]
 
tss$chromosome<-as.character(tss$chromosome)
w<-which(tss$chromosome=="chrX")
tss[w,"chromosome"]<-"chr23"
 
w<-which(tss$chromosome=="chrY")
tss[w,"chromosome"]<-"chr24"

w<-which(tss$V7=="-")
tss$start.original<-tss$start
tss$end.original<-tss$end
 
tss$start[w]<-tss$end[w]
tss$end[w]<-tss$start.original[w]
 
tss$end<-tss$start
 
f.all<-data.frame()
 
#"chr7:5336538:SNP"
##
for (i in 1:23){
cat(i,"\n")
cl_chr<-cl[which(cl$snpchr%in%paste0("chr",i)),]
a <- read_tsv(paste0("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.filtered.",i,".ld.formatted.gz"))
a<-data.frame(a)
w<-which(names(a)%in%"X11")
if(length(w)>0){a<-a[,-w]}

#a$BP_A<-as.character(a$BP_A)
#a$BP_B<-as.character(a$BP_B)
#a$CHR_A<-as.character(a$CHR_A)
#a$CHR_B<-as.character(a$CHR_B)

agg1<-a %>%  
group_by(SNP_A) %>%  
dplyr::summarise (bpa = dplyr::first(BP_A),start_bp = min(BP_B),  
              stop_bp = max(BP_B),nproxies = n(),chr=dplyr::first(CHR_A))
 
agg2<-a %>%  
   group_by(SNP_B) %>%  
   dplyr::summarize (bpb = dplyr::first(BP_B),start_bp = min(BP_A),  
              stop_bp = max(BP_A),nproxies = n(),chr=dplyr::first(CHR_B))
 
b <- full_join(agg1, agg2, c("SNP_A" = "SNP_B"))
o<-order(as.numeric(b$bpa))
b<-b[o,]
 
#b<-merge(agg1,agg2,by.x="SNP_A",by.y="SNP_B",all=T)
 
 
 #b$final_min <- b$start_bp.x
 #b$final_max <- b$stop_bp.y
 #b$final_min[is.na(b$final_min)] <- b$bpa[is.na(b$final_min)]
 #b$final_max[is.na(b$final_max)] <- b$bpb[is.na(b$final_max)]
 #b<-data.frame(b)
 
b2<-data.frame(b$nproxies.x,b$nproxies.y)
b$nproxies<-apply(b2,1,function(x) y<-sum(x,na.rm=T))
 
b2<-data.frame(b$start_bp.x,b$start_bp.y,b$bpa)
b$final_min<-apply(b2,1,function(x) y<-min(x,na.rm=T))
 
b2<-data.frame(b$stop_bp.x,b$stop_bp.y,b$bpa)
b$final_max<-apply(b2,1,function(x) y<-max(x,na.rm=T))
 
 
b$CHR<-b$chr.x
b$CHR[is.na(b$CHR)]<-b$chr.y[is.na(b$CHR)]
b$CHR<-paste0("chr",b$CHR)
 
table(is.na(b$final_min))
table(is.na(b$final_max))
table(b$final_min <= b$final_max)
table(b$final_min <= b$final_max)
 
 
b <- dplyr::select(b, CHR=CHR,SNP=bpa, min=final_min, max=final_max,nproxies=nproxies)
summary(b$max - b$min)



#add MAF
f <- read_tsv(paste0("/panfs/panasas01/shared-godmc/1kg_reference_ph3/frq/eur.filtered.",i,".formatted.gz"))
spl<-strsplit(f$SNP,split=":")
spl<-do.call("rbind",spl)
colnames(spl)<-c("snpchr","snppos","snptype")
f<-data.frame(spl,f)
m<-match(f$snppos,b$SNP)
f<-data.frame(snpchr=f[,c("snpchr")],b[m,-1:-2],f[,c("snppos","MAF","snptype","SNP")])
f$snppos<-as.numeric(as.character(f$snppos))
f$min<-as.numeric(as.character(f$min))
f$max<-as.numeric(as.character(f$max))
 
f$min[is.na(f$min)]<-f$snppos[is.na(f$min)]
f$max[is.na(f$max)]<-f$snppos[is.na(f$max)]
f$nproxies[is.na(f$nproxies)]<-0


#add distance to CpG
#load("../03_clumping_16/cpg_pos.rdata")
#cpgpos<-cpgpos[which(cpgpos$cpgchr==paste0("chr",i)),]
 
 
#cpgpos<-y[which(y$chromosome==paste0("chr",i)),]
if(i==23){
j<-"X"
cpgpos<-y[which(y$chromosome==paste0("chr",j)),]
cpgpos$chromosome<-"chr23"}
 
if(i<23){cpgpos<-y[which(y$chromosome==paste0("chr",i)),]}
 
cpgpos$position<-as.numeric(as.character(cpgpos$position))
snppos<-as.numeric(as.character(f$snppos))


#distance to cpg
 
#Matt's function is really fast, don't use this.
#closestcpg<-ddply(f, .(snppos), summarise,
#       closestcpgposmin = snppos-min(abs(snppos-cpgpos$position),na.rm=T),closestcpgposmax = snppos+min(abs(snppos-cpgpos$position),na.rm=T) )
 
#m1<-match(closestcpg[,2],cpgpos$position)
#m2<-match(closestcpg[,3],cpgpos$position)
#cpgpos1<-data.frame(cpgpos[m1,"name"],cpgpos[m2,"name"])
 
#cpgpos2<-as.character(cpgpos1[,1])
#cpgpos2[!is.na(cpgpos1[,2])]<-as.character(cpgpos1[!is.na(cpgpos1[,2]),2])
#m<-match(cpgpos2,cpgpos$name)
#dist<-abs(cpgpos[m,c("position")]-f$snppos)
#f<-data.frame(f,closest450kcpg=cpgpos[m,c("name")],closest450kcpgpos=cpgpos[m,c("position")],closest450kdistance=dist) #16050654      72612
 
#distance to CpG
 
f$chromosome<-f$snpchr
f$position<-f$snppos
cpgpos$start<-cpgpos$position
cpgpos$end<-cpgpos$position
 
 
res<-genomic.nearest(positions=f,intervals=cpgpos)
res<-cpgpos[res,c("name","position")]
cpgdist<-res$position-f$snppos
 
f<-data.frame(f,closest450kcpg=res$name,closest450kcpgpos=res$position,closest450kdistance=cpgdist) #16050654      72612
 
#distance to TSS
res<-genomic.nearest(positions=f,intervals=tss)
res<-tss[res,]
spl<-do.call("rbind",strsplit(as.character(res$V9),split=";"))
tssdist<-res$start-f$snppos
genename<-gsub("gene_name ","",spl[,4])
f<-data.frame(f,tssdist,closesttss=res$start,closestgene=genename)

###

#add snp_cis 
snp_cis<-read_tsv(paste0("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr",i,".cistrans3.txt.gz"))

snp_cis<-unique(snp_cis[,c("snp","snp_cis","min_pval","max_abs_Effect","max_abs_SE","TotalSampleSize","Freq1","HetISq")])
names(snp_cis)<-c("snp_tested","snp_cis","min_pval","max_abs_Effect","max_abs_SE","TotalSampleSize","Freq1","HetISq")

m<-match(f$SNP,snp_cis$snp_tested)
f<-data.frame(f,snp_cis[m,])

f.all<-rbind(f.all,f)
 
}

#add clumped mqtl

f.all$mqtl_clumped<-f.all$snp_tested
w<-which(!is.na(f.all$mqtl_clumped))
f.all$mqtl_clumped[w]<-"FALSE"

w<-which(f.all$SNP%in%cl.snp$mqtl_clumped)
f.all$mqtl_clumped[w]<-"TRUE"

print(table(f.all$mqtl_clumped))

save(f.all,file="../results/enrichments/snpcontrolsets_se.rdata")





