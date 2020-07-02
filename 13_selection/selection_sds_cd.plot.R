library(ggplot2)
library(dplyr)
library(data.table)

contingency<-function (af, prop, odds_ratio, eps = 1e-15) 
{
    a <- odds_ratio - 1
    b <- (af + prop) * (1 - odds_ratio) - 1
    c_ <- odds_ratio * af * prop
    if (abs(a) < eps) {
        z <- -c_/b
    }
    else {
        d <- b^2 - 4 * a * c_
        if (d < eps * eps) {
            s <- 0
        }
        else {
            s <- c(-1, 1)
        }
        z <- (-b + s * sqrt(max(0, d)))/(2 * a)
    }
    y <- vapply(z, function(a) zapsmall(matrix(c(a, prop - a, 
        af - a, 1 + a - af - prop), 2, 2)), matrix(0, 2, 2))
    i <- apply(y, 3, function(u) all(u >= 0))
    return(y[, , i])
}

get_population_allele_frequency<-function (af, prop, odds_ratio, prevalence) 
{
    co <- contingency(af, prop, odds_ratio)
    af_controls <- co[1, 2]/(co[1, 2] + co[2, 2])
    af_cases <- co[1, 1]/(co[1, 1] + co[2, 1])
    af <- af_controls * (1 - prevalence) + af_cases * prevalence
    return(af)
}


get_r_from_lor <- function(lor, af, ncase, ncontrol, prevalence, model="logit")
{
        stopifnot(length(lor) == length(af))
        stopifnot(length(ncase) == 1 | length(ncase) == length(lor))
        stopifnot(length(ncontrol) == 1 | length(ncontrol) == length(lor))
        stopifnot(length(prevalence) == 1 | length(prevalence) == length(lor))
        if(length(prevalence) == 1 & length(lor) != 1)
        {
                prevalence <- rep(prevalence, length(lor))
        }
        if(length(ncase) == 1 & length(lor) != 1)
        {
                ncase <- rep(ncase, length(lor))
        }
        if(length(ncontrol) == 1 & length(lor) != 1)
        {
                ncontrol <- rep(ncontrol, length(lor))
        }

        nsnp <- length(lor)
        r <- array(NA, nsnp)
        for(i in 1:nsnp)
        {
                if(model == "logit")
                {
                        ve <- pi^2/3
                } else if(model == "probit") {
                        ve <- 1
                } else {
                        stop("Model must be probit or logit")
                }
                popaf <- get_population_allele_frequency(af[i], ncase[i] / (ncase[i] + ncontrol[i]), exp(lor[i]), prevalence[i])
                vg <- lor[i]^2 * popaf * (1-popaf)
                r[i] <- sqrt(vg / (vg + ve) / 0.58)
        }
        return(r)
}

prevalence=0.005
ncase=12882
ncontrol=21770

#r2<-read.table("./cd_sds/vars.IBD_CD_mqtl_6cols",sep=" ",he=T)
r2<-read.table("./cd_sds/vars.IBD_CD_Liu_2015",sep=" ",he=T)
r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./cd_sds/vars.IBD_CD_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,8]),":SNP")

sds<-paste0("chr",r3[,8],":SNP")


##ES = 2β^2f(1 − f)
#f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
#f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))

#p1<-ggplot(f.all, aes(es, colour=cd_mqtl)) +
#geom_density() +
#labs(x="Genetic Variance")
#ggsave(p1,file="Mqtl_GeneticVariance.pdf")

#how much of the cd variance
h<-fread("./cd_sds/EUR.IBD.gwas_info03_filtered.assoc",he=T)
w1<-which(h$A1=="D")
w2<-which(h$A2=="D")
w<-unique(c(w1,w2))
h$SNP<-"SNP"
h$SNP[w]<-"INDEL"
id<-paste0("chr",h$CHR,":",h$BP,":",h$SNP)
h<-data.frame(h,id)

#bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
#m<-match(h$SNP,bim$V2)
#h<-data.frame(SNP=paste0("chr",bim[m,1],":",bim[m,4],":SNP"),h)

#h$b_sq<-h$OR^2
h$MAF<-h$FRQ_U_21770
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
#h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))
h$logodds<-log(h$OR)
h$es<-get_r_from_lor (lor=h$logodds, af=h$MAF, ncase=ncase, ncontrol=ncontrol, prevalence=prevalence, model="logit")
save(h,file="ibd.Robj")
#r<-read.table("./cd_sds/cd_snps_chrpos.txt",he=F,sep="\t")

cd<-fread("./cd_sds/EUR.CD.gwas_info03_filtered.assoc")
r<-read.table("./cd_sds/cd_liu2016.txt",sep="\t",he=T)
r$id<-paste0("chr",r$CHR,":",r$BP,":SNP")
w<-which(r$A1=="D"|r$A2=="D")
r$id[w]<-paste0("chr",r$CHR[w],":",r$BP[w],":INDEL")


m<-match(r$SNP,cd$SNP)
cd2<-cd[m,]

cd2<-na.omit(cd2[which(is.na(r$GWAS.P)),])
add.SNP<-paste0(cd2$CHR,":",cd2$BP,":SNP")
r<-r[which(r$GWAS.P<5e-8),]
#r<-read.csv("./cd_sds/studies_GCST003043-associations-2020-02-14.csv")
#r<-read.csv("./cd_sds/publications_26192919-associations-2020-02-14.csv")
#r<-r[which(r$Reported.trait=="Crohn's disease"),]
#spl<-do.call("rbind",strsplit(as.character(r[,1]),split="-"))

#write.table(spl[,1],"./cd_sds/IBD.txt",sep="\t",quote=F,row.names=F,col.names=F)

#r<-data.frame(rsid=spl[,1],r)
#r_pos<-read.csv("./cd_sds/IBD_pos.txt",sep="\t",he=T)
#m<-match(r$rsid,r_pos$name)
#r<-data.frame(r,r_pos[m,])
#r$id<-paste0(r$chrom,":",r$chromEnd,":SNP")
#spl2<-do.call("rbind",strsplit(as.character(r$P.value),split=" x "))
#spl2[,2]<-gsub("10-","e-",spl2[,2])
#r$Pval<-paste0(spl2[,1],spl2[,2])
#r$Pval<-as.numeric(r$Pval)
#r<-r[which(r$Pval<5e-8),]

h$cd_mqtl<-"no_mqtl"
w<-which(h$SNP%in%c(r$id,add.SNP))
h$cd_mqtl[w]<-"cd SNPs (n=31)"
h_all<-h[w,]

w<-which(h$SNP%in%mqtl)
h$cd_mqtl[w]<-"cd mqtl (n=60)"
h_all2<-h[w,]

w<-which(h$SNP%in%sds)
h$cd_mqtl[w]<-"cd mqtl sds (n=3)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)
#h_all<-rbind(h_all,h_all3)

table(h_all$cd_mqtl)

#   cd mqtl (n=60) cd mqtl sds (n=3)    cd SNPs (n=31) 
#               60                 3                31 

#   cd mqtl (n=60) cd mqtl sds (n=3)    cd SNPs (n=31) 
#               91                 2               199 


p1<-ggplot(h_all, aes(es, colour=cd_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="cd_GeneticVariance.pdf")


h_all%>%group_by(cd_mqtl)%>%summarise(es=mean(es,na.rm=T))
# A tibble: 3 x 2
#            cd_mqtl         es
#              <chr>      <dbl>
#1    cd mqtl (n=60) 0.06206499
#2 cd mqtl sds (n=3) 0.06075588
#3    cd SNPs (n=31) 0.03340113

h_all$cd_mqtl <- factor(h_all$cd_mqtl, levels = c("cd SNPs (n=31)","cd mqtl (n=60)", "cd mqtl sds (n=3)"))

#no mapping
load("/newshared/godmc/database_files/snps.rdata")
h_all[which(h_all$SNP%in%out_df2$name==F),c("SNP","cd_mqtl")]
#                      SNP        cd_mqtl
#192126  chr10:35542343:SNP cd mqtl (n=60)
#3708494 chr16:50756540:SNP cd mqtl (n=60)
#8699915  chr6:32767249:SNP cd mqtl (n=60)

length(mqtl) #108
h_all[which(h_all$SNP%in%mqtl==F),c("SNP","cd_mqtl")]
length(unique(h_all[which(h_all$SNP%in%mqtl),c("SNP")])) #91
unique(h_all[which(h_all$SNP%in%mqtl),c("SNP")])

length(sds)
h_all[which(h_all$SNP%in%sds==F),c("SNP","cd_mqtl")]
length(unique(h_all[which(h_all$SNP%in%sds),c("SNP")])) #2
unique(h_all[which(h_all$SNP%in%sds),c("SNP")])

p1<-ggplot(h_all, aes(y=es, x=cd_mqtl)) +
  geom_boxplot() +
  labs(y="Genetic Variance",x="mQTL category")
ggsave(p1,file="CD_GeneticVariance_boxplot.pdf")

save(h_all,file="./cd_sds/cd_plot.Robj")

