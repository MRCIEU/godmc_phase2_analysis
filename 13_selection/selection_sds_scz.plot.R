library(ggplot2)
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
spl<-strsplit(as.character(h$SNP),split=":")
spl<-do.call("rbind",spl)

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

load("../results/enrichments/snpcontrolsets_selection.rdata")
w<-which(is.na(f.all$snp_cis))
f.all$Category<-as.character(f.all$snp_cis)
f.all$Category[w]<-"no_mqtl"

f.all$Category<-gsub("TRUE","cisonly",f.all$Category)
f.all$Category<-gsub("FALSE","transonly",f.all$Category)
f.all$Category<-gsub("ambivalent","cis+trans",f.all$Category)

f.all$min_log10pval<-f.all$min_pval
w0<-which(f.all$min_pval==0)
mx<-min(f.all$min_pval[-w0],na.rm=T)
f.all$min_log10pval[w0]<-mx
f.all$min_log10pval<--log10(as.numeric(f.all$min_log10pval))

w<-which(f.all$mqtl_clumped=="TRUE")
f.all2<-f.all[w,]

r2<-read.table("./scz_sds/vars.SCZ_Pardinas_2018_mqtl_7cols",sep=" ",he=T)
r2[,8]<-gsub("\\([^\\)]+\\)","",as.character(r2[,5])) #260 #40

r3<-read.table("./scz_sds/vars.SCZ_Pardinas_2018_sds",sep=" ",he=T)
r3[,8]<-gsub("\\([^\\)]+\\)","",as.character(r3[,5]))

mqtl<-paste0("chr",unique(r2[,8]),":SNP")

f.all$scz_mqtl<-"no mqtl"

w<-which(f.all$mqtl_clumped=="TRUE")
f.all$scz_mqtl[w]<-"clumped mqtl"

w<-which(f.all$SNP%in%mqtl)
f.all$scz_mqtl[w]<-"scz mqtl"

sds<-paste0("chr",r3[,8],":SNP")
w<-which(f.all$SNP%in%sds)
w<-which(sds%in%f.all$SNP==F)
if(length(w)>0){
sds[w]<-gsub(":SNP",":INDEL",sds[w])}

w<-which(f.all$SNP%in%sds)
f.all$scz_mqtl[w]<-"scz mqtl sds"

##ES = 2β^2f(1 − f)
#f.all$max_abs_Effect_sq<-f.all$max_abs_Effect^2
#f.all$es<-(2*f.all$max_abs_Effect_sq)*(f.all$MAF*(1-f.all$MAF))

p1<-ggplot(f.all, aes(es, colour=scz_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="Mqtl_GeneticVariance.pdf")


r<-read.table("./scz_sds/scz_gwas.txt",he=T)
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig",sep="\t",he=F)
m<-match(r$rsid,bim$V2)
r<-data.frame(r,pos=bim$V4[m])
r$SNP<-paste0("chr",r$chr,":",r$pos,":",r$INDEL)

w<-which(is.na(r$pos))#17
r<-r[-w,]

h$scz_mqtl<-"No mqtl"
w<-which(h$ID%in%r$SNP)
h$scz_mqtl[w]<-"scz SNPs (n=159)"
h_all<-h[w,]

w<-which(h$ID%in%mqtl)
h$scz_mqtl[w]<-"scz mQTL (n=133)"
h_all2<-h[w,]

w<-which(h$ID%in%sds)
h$scz_mqtl[w]<-"scz mQTL SDS (n=9)"
h_all3<-h[w,]

h_all<-rbind(h_all,h_all2,h_all3)
#h_all<-rbind(h_all,h_all3)

table(h_all$scz_mqtl)

#scz mQTL (n=133) scz mQTL SDS (n=9)   scz SNPs (n=159)
#133                  9                159

p1<-ggplot(h_all, aes(es, colour=scz_mqtl)) +
geom_density() +
labs(x="Genetic Variance")
ggsave(p1,file="scz_GeneticVariance.pdf")

library(dplyr)
h_all%>%group_by(scz_mqtl)%>%summarise(es=mean(es,na.rm=T))
# A tibble: 3 x 2
#            cd_mqtl         es
#              <chr>      <dbl>
#1    cd mqtl (n=60) 0.06206499
#2 cd mqtl sds (n=3) 0.06075588
#3    cd SNPs (n=31) 0.03340113

h_all$scz_mqtl <- factor(h_all$scz_mqtl, levels = c("scz SNPs (n=159)","scz mQTL (n=133)", "scz mQTL SDS (n=9)"))
p1<-ggplot(h_all, aes(y=es, x=scz_mqtl)) +
  geom_boxplot() +
  labs(y="Genetic Variance",x="mQTL category")
ggsave(p1,file="scz_GeneticVariance_boxplot.pdf")
save(h_all,file="./scz_sds/scz_plot.Robj")


save(h_all,file="./scz_sds/scz_plot.Robj")

