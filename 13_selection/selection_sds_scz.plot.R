library(ggplot2)
library(dplyr)

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

prevalence=0.01
ncase=40675
ncontrol=64643

h<-read.table("./scz_sds/clozuk_pgc2.meta.sumstats.txt.gz",he=T)

h<-data.frame(ID=paste0("chr",h$CHR,":",h$BP,":","SNP"),h)
a1<-(nchar(as.character(h$A1)))
a2<-(nchar(as.character(h$A2)))
w1<-which(a1>1)
w2<-which(a2>1)
length(unique(c(w1,w2)))
w<-unique(c(w1,w2))
id<-paste0("chr",h$CHR[w],":",h$BP[w],":","INDEL")
h$ID<-as.character(h$ID)
h$ID[w]<-id

h$MAF<-h$Freq.A1
w<-which(h$MAF>0.5)
h$MAF[w]<-1-h$MAF[w]
#h$es<-(2*h$b_sq)*(h$MAF*(1-h$MAF))
h$logodds<-log(h$OR)
h$es<-get_r_from_lor(lor=h$logodds, af=h$MAF, ncase=ncase, ncontrol=ncontrol, prevalence=prevalence, model="logit")
save(h,file="./scz_sds/clozuk_pgc2.Robj")


#r2<-read.table("./scz_sds/vars.SCZ_Pardinas_2018_mqtl_7cols",sep=" ",he=T)
r2<-read.table("./scz_sds/vars.SCZ_Pardinas_2018",sep=" ",he=T)

r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./scz_sds/vars.SCZ_Pardinas_2018_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,8]),":SNP")
sds<-paste0("chr",r3[,8],":SNP")


##ES = 2β^2f(1 − f)
#f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
#f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))


r<-read.table("./scz_sds/scz_gwas.txt",he=T) # no match in dbSNP for rs72342102
r_pos<-read.table("./scz_sds/scz_pos.txt",he=T) 

m<-match(r$rsid,r_pos$name)
r<-data.frame(r,r_pos[m,])
r$id<-paste0(r$chrom,":",r$chromEnd,":",r$INDEL)

w<-which(is.na(r$chromEnd)) #3
r<-r[-w,]

h$scz_mqtl<-"No mqtl"
w<-which(h$ID%in%r$id)
h$scz_mqtl[w]<-"scz SNPs (n=159)"
h_all<-h[w,]

w<-which(h$ID%in%mqtl)
#h$scz_mqtl[w]<-"scz mQTL (n=133)"
h$scz_mqtl[w]<-"scz mQTL (n=226)"
h_all2<-h[w,]

w<-which(h$ID%in%sds)
h$scz_mqtl[w]<-"scz mQTL SDS (n=9)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)
#h_all<-rbind(h_all,h_all3)

table(h_all$scz_mqtl)

#scz mQTL (n=133) scz mQTL SDS (n=9)   scz SNPs (n=159)
#133                  9                159

table(h_all$scz_mqtl)

#  scz mQTL (n=226) scz mQTL SDS (n=9)   scz SNPs (n=159) 
#               226                  9                165

p1<-ggplot(h_all, aes(es, colour=scz_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="scz_GeneticVariance.pdf")


h_all%>%group_by(scz_mqtl)%>%summarise(es=mean(es,na.rm=T))
# A tibble: 3 x 2
#  scz_mqtl               es
#  <chr>               <dbl>
#1 scz mQTL (n=226)   0.0245
#2 scz mQTL SDS (n=9) 0.0298
#3 scz SNPs (n=159)   0.0225

h_all$scz_mqtl <- factor(h_all$scz_mqtl, levels = c("scz SNPs (n=159)","scz mQTL (n=226)", "scz mQTL SDS (n=9)"))
p1<-ggplot(h_all, aes(y=es, x=scz_mqtl)) +
  geom_boxplot() +
  labs(y="Genetic Variance",x="mQTL category")
ggsave(p1,file="scz_GeneticVariance_boxplot.pdf")
save(h_all,file="./scz_sds/scz_plot.Robj")

#no mapping
load("/newshared/godmc/database_files/snps.rdata")
h_all[which(h_all$ID%in%out_df2$name==F),c("SNP","scz_mqtl")]
#                      SNP        cd_mqtl
#192126  chr10:35542343:SNP cd mqtl (n=60)
#3708494 chr16:50756540:SNP cd mqtl (n=60)
#8699915  chr6:32767249:SNP cd mqtl (n=60)

length(mqtl) #263
h_all[which(h_all$ID%in%mqtl==F),c("SNP","scz_mqtl")]
length(unique(h_all[which(h_all$ID%in%mqtl),c("SNP")])) #226
unique(h_all[which(h_all$ID%in%mqtl),c("SNP")])

length(sds) #9
h_all[which(h_all$ID%in%sds==F),c("SNP","scz_mqtl")]
length(unique(h_all[which(h_all$ID%in%sds),c("SNP")])) #9
unique(h_all[which(h_all$ID%in%sds),c("SNP")])


