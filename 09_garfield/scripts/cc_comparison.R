#It could be one of these three models
 
#G -> cpg -> cell counts
#G -> cell counts -> cpg
#Cell counts <- G -> cpg
 
#I think we could look at a couple of things using the mQTLs
 
#abs(beta_cpg) ~ abs(beta_cellcount)
#proportion for which r^2(g->cpg) > r^2(g->cellcount)

#R2 = F / (F + n – 2)
 #Where F = b^2 / se^2 

#For model (A) would expect the result from (2) to have large proportion of r2 being higher for cpg than cellcounts

library(data.table)
library(ggplot2)

arguments<-commandArgs(T)
trait<-as.character(arguments[1])

#trait<-as.character("27863252-GCST004599-EFO_0004584")
spl<-strsplit(trait,split="-")
trait_id<-spl[[1]][2]

#path="/projects/MRC-IEU/research/projects/ieu2/p5/021/working/data/Astle2016/"
#meta<-read.table(paste0(path,"meta_data_forMRBase_Astle2016.txt"),he=T,sep="\t")

meta<-read.table("meta_data_forMRBase_Astle2016.txt",he=T,sep="\t")
id<-as.character(meta[which(meta$Study.accession%in%trait_id),c("Reported.trait")])
id<-gsub(" ","_",id)

cc.n<-as.character(meta[which(meta$Study.accession%in%trait_id),c("samplesize")])

#l<-list.files(path=path,pattern=paste0("27863252-",id))

#cc<-read.table(paste0(path,l),he=T)
#cc<-fread(paste0(path,l),he=T,sep="\t")
#cc<-fread("gzip -dc GCST004631.txt",he=T,sep=" ")

cc<-fread(trait,he=T,sep=" ")
cc$id<-paste0("chr",cc$chromosome,":",cc$base_pair_location,":SNP")
w1<-which(nchar(cc$effect_allele)>1)
w2<-which(nchar(cc$other_allele)>1)
w<-unique(c(w1,w2))
cc$id[w]<-gsub("SNP","INDEL",cc$id[w])

load("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped<- subset(clumped, pval < 1e-8 & cis | pval < 1e-14 & !cis)

m<-match(clumped$snp,cc$id)
cc2<-cc[m,]

length(which(is.na(cc2$id)))
w<-which(is.na(cc2$id))
cc2<-cc2[-w]
clumped<-clumped[-w,]

cc2$a1<-cc2$effect_allele
cc2$a2<-cc2$other_allele

w1<-which(nchar(cc2$a2)>nchar(cc2$a1))
w2<-which(nchar(cc2$a2)<nchar(cc2$a1))

cc2$a2[w1]<-"I"
cc2$a1[w1]<-"D"

cc2$a2[w2]<-"D"
cc2$a1[w2]<-"I"

clumped$Allele1 <- toupper(clumped$Allele1)
clumped$Allele2 <- toupper(clumped$Allele2)
table(clumped$Allele1 == cc2$a1 | clumped$Allele2 == cc2$a2)

w <- which(clumped$Allele1 == cc2$a2 & clumped$Allele2==cc2$a1)
cc2$beta_a1<- cc2$beta
cc2$beta_a1[w]<- -(cc2$beta[w])
cc2$a1_freq<-cc2$effect_allele_frequency
cc2$a1_freq[w]<-1-cc2$effect_allele_frequency[w]

cc2$a1_2<-cc2$a1
cc2$a2_2<-cc2$a2

cc2$a1[w]<-cc2$a2_2[w]
cc2$a2[w]<-cc2$a1_2[w]

table(clumped$Allele1 == cc2$a1 | clumped$Allele2 == cc2$a2)
# TRUE 
#267382

#w<-which(clumped$Allele1 == cc2$a1 | clumped$Allele2 == cc2$a2)
#head(clumped[-w,c("snp","Allele1","Allele2","Freq1")])
#head(cc2[-w,c("id","a1","a2","effect_allele_frequency")])
df<-data.frame(clumped,cc2)
w<-which(df$p_value<5e-8)
df$GS<-"P>5e-8"
df$GS[w]<-"P<5e-8"

df2<-df[which(df$p_value<0.2),]

p1<-ggplot(df, aes(x=Effect, y=beta_a1)) +
geom_point(aes(colour=GS), cex=0.5) +
geom_point(data=subset(df, GS=="P<5e-8"), aes(x = Effect, y = beta_a1, color = GS), cex=0.5) +

labs(x="GoDMC effect size", y=paste0(id," effect size"), colour="GWA Pval") +
theme(legend.position="right") +
scale_fill_brewer(type="qual") +
guides(alpha=FALSE)
ggsave(p1,file=paste0("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/09_garfield/images/",id,".png"), width=9, height=4)

##
df$F <- df$Effect^2 / df$StdErr^2
df$r2 <- df$F / (df$F + df$TotalSampleSize - 2)

df$beta_sq<-df$beta_a1^2
df$se_sq<-df$standard_error^2
df$F_cc<-df$beta_a1^2/df$standard_error^2
df$r2_cc<-df$F_cc/(df$F_cc+as.numeric(cc.n)-2)
df$p_value <- as.numeric(df$p_value)

df1 <- subset(df, p_value < 0.05)

df2 <- subset(df, p_value < (0.05/nrow(df)))
cat("Number of overlaps: ", nrow(df),"\n")
cat(length(which(df$r2 > df$r2_cc))/nrow(df),"\n")
cat("Number of cc with p < 0.05: ", nrow(df1),"\n")
cat(length(which(df1$r2 > df1$r2_cc))/nrow(df1),"\n")
cat("Number of cc with p < 0.05/ntest: ", nrow(df2),"\n")
cat(length(which(df2$r2 > df2$r2_cc))/nrow(df2),"\n")

df3<-data.frame(id=id,r2_godmc_larger_than_cc=length(which(df2$r2 > df2$r2_cc)),no_assoc=nrow(df),prop=length(which(df2$r2 > df2$r2_cc))/nrow(df))
write.table(df3,paste0("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/09_garfield/cellcounts/",id,".txt"),sep="\t",quote=F,row.names=F,col.names=T)
#save(df,file=paste0("/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/09_garfield/cellcounts/",id,".Robj"))


