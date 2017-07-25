library(tidyverse)
library(dplyr)
library(meffil)

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
retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]
#420509
 
nrow(cl)

#316245
cl<-cl[which(cl$cpg%in%retaincpg),]
nrow(cl)
#[1] 249008

 
length(unique(cl$snp))
#[1] 186395

y<-meffil.get.features("450k")


#get TSS data
#distance to TSS

tss<-read.table("gencode.v26.annotation.gtf.gz",skip=5,sep="\t")
names(tss)[1:5]<-c("chromosome","source","type","start","end")

table(tss$type)



#       CDS           exon           gene Selenocysteine    start_codon 
#        710618        1194547          58219            119          83245 
#    stop_codon     transcript            UTR 
#        74994         199324         283420 
 
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
 
agg1<-a %>%  
group_by(SNP_A) %>%  
dplyr::summarize (bpa = first(BP_A),start_bp = min(BP_B),  
              stop_bp = max(BP_B),nproxies = n(),chr=first(CHR_A))
 
agg2<-a %>%  
   group_by(SNP_B) %>%  
   dplyr::summarize (bpb = first(BP_B),start_bp = min(BP_A),  
              stop_bp = max(BP_A),nproxies = n(),chr=first(CHR_B))
 
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
cpgdist<-abs(res$position-f$snppos)
 
f<-data.frame(f,closest450kcpg=res$name,closest450kcpgpos=res$position,closest450kdistance=cpgdist) #16050654      72612
 
#distance to TSS
res<-genomic.nearest(positions=f,intervals=tss)
res<-tss[res,]
tssdist<-abs(res$start-f$snppos)
f<-data.frame(f,tssdist,closesttss=res$start)

###
cl_chr.cis<-cl_chr[which(cl_chr$cis==T),]
cl_chr.trans<-cl_chr[which(cl_chr$cis==F),]

m<-which(f$SNP%in%cl_chr.cis$snp)
f$cismQTL<-"FALSE"
f$cismQTL[m]<-"TRUE"
 
m<-which(f$SNP%in%cl_chr.trans$snp)
f$transmQTL<-"FALSE"
f$transmQTL[m]<-"TRUE"
 
m<-which(f$SNP%in%unique(cl_chr$snp))
f$mQTL<-"FALSE"
f$mQTL[m]<-"TRUE"


##make categories

mafcat<-cut(f$MAF,breaks=seq(0,0.5,0.05))
 
cpgdistcat<-cut(f$closest450kdistance,breaks=seq(0,10000000,10000))
w<-which(f$closest450kdistance==0)
cpgdistcat[w]<-levels(cpgdistcat)[1]
 
proxycat<-cut(f$nproxies,breaks=seq(0,100000,10))
w<-which(f$nproxies==0)
proxycat[w]<-levels(proxycat)[1]
 
tssdistcat500<-cut(f$tssdist,breaks=seq(0,10000000,500))
w<-which(f$tssdist==0)
tssdistcat500[w]<-levels(tssdistcat500)[1]
 
tssdistcat2000<-cut(f$tssdist,breaks=seq(0,10000000,2000))
w<-which(f$tssdist==0)
tssdistcat2000[w]<-levels(tssdistcat2000)[1]
 
tssdistcat5000<-cut(f$tssdist,breaks=seq(0,10000000,5000))
w<-which(f$tssdist==0)
tssdistcat5000[w]<-levels(tssdistcat5000)[1]
 
tssdistcat10000<-cut(f$tssdist,breaks=seq(0,10000000,10000))
w<-which(f$tssdist==0)
tssdistcat10000[w]<-levels(tssdistcat10000)[1]
 
tssdistcat20000<-cut(f$tssdist,breaks=seq(0,10000000,20000))
w<-which(f$tssdist==0)
tssdistcat20000[w]<-levels(tssdistcat20000)[1]
 
tssdistcat100000<-cut(f$tssdist,breaks=seq(0,10000000,100000))
w<-which(f$tssdist==0)
tssdistcat100000[w]<-levels(tssdistcat100000)[1]
 
f<-data.frame(f,mafcat,cpgdistcat,proxycat,tssdistcat500bp=tssdistcat500,tssdistcat2kb=tssdistcat2000,tssdistcat5kb=tssdistcat5000,tssdistcat10kb=tssdistcat10000,tssdistcat20kb=tssdistcat20000,tssdistcat100kb=tssdistcat100000)

f$groups<-paste(f$mafcat,f$cpgdistcat,f$nproxies,f$tssdistcat500bp)
#f$groups2<-paste(f$mafcat,f$cpgdistcat,f$proxycat,f$tssdistcat500bp)
#f$groups3<-paste(f$mafcat,f$cpgdistcat,f$proxycat,f$tssdistcat2kb)
#f$groups4<-paste(f$mafcat,f$cpgdistcat,f$proxycat,f$tssdistcat5kb)
#f$groups5<-paste(f$mafcat,f$cpgdistcat,f$proxycat,f$tssdistcat10kb)
#f$groups6<-paste(f$mafcat,f$cpgdistcat,f$proxycat,f$tssdistcat20kb)
 
 
f.all<-rbind(f.all,f)
save(f.all,file="../results/enrichments/snpcontrolsets.rdata")
 
}


