regional_plot<-function(diab_coloc=diab_coloc,snpid=snpid,cpgid=cpgid,trait=trait,gwa=gwa,mQTL=mQTL)
{
snppos<-diab_coloc[diab_coloc$V2==snpid,"V4"]/1e6
snppos_bp<-diab_coloc[diab_coloc$V2==snpid,"V4"]
chr<-diab_coloc[diab_coloc$V2==snpid,"V1"]
start_bp<-(snppos_bp-1e6)
stop_bp<-(snppos_bp+1e6)
start<-start_bp/1e6
stop<-stop_bp/1e6
GWAS<-gwa[which(gwa$V1==chr & gwa$V4>start_bp& gwa$V4<stop_bp),]
GWAS$BP<-GWAS$V4/1e6
GWAS$log10P<--log10(GWAS$P)
gwa_p<-GWAS[GWAS$SNP==snpid,"log10P"]

#mQTL<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr11.cistrans.txt.gz")
mQTL$BP<-mQTL$snppos/1e6
mQTL$log10P<--log10(mQTL$pval)
w<-which(mQTL$pval==0)
mQTL$log10P[w]<--log10(min(mQTL$pval[-w]))
mqtl_p<-mQTL[mQTL$snppos==snppos_bp,"log10P"]

pdf(paste0(trait,"_",cpgid,"_",snpid,".pdf"),height=6,width=12)
par(mar = c(5, 5, 3, 5))
with(GWAS, plot(BP, log10P, pch=16, cex=1.2, cex.lab=1.2, col="grey", ylim=c(0,max(GWAS$log10P)+1), xlim=c(start,stop),
                ylab=expression(-log[10](italic(p))[trait]),
                xlab="Genomic position (Mb)"))
par(new = T)
points(snppos,gwa_p,pch=24,col="black",xlim=c(start,stop),ylim=c(0,max(GWAS$log10P)+1))
text(snppos-0.1,gwa_p+0.3,snpid,col="black")
par(new = T)
with(mQTL, plot(BP, log10P, pch=18, axes=F, xlab=NA, ylab=NA, cex=1.2, cex.lab=1.2, xlim=c(start,stop), col="lightgreen"))
par(new = T)
points(snppos,mqtl_p,pch=24,col="black",xlim=c(start,stop),ylim=c(0,max(GWAS$log10P)+1))
text(snppos-0.1,mqtl_p-15,snpid,col="darkgreen")
axis(side = 4)
mtext(side = 4, line = 3, expression(-log[10](italic(p))[Methylation]))
legend(start,max(mQTL$log10P)-10,
       legend=c(trait, cpgid),
       lty=c(0,0,0), pch=c(16, 18,18), col=c("grey","lightgreen"))
dev.off()
}

snpid="rs2490852"
cpgid="cg14451791"
trait="Fasting_insulin"

load("../data/diab_coloc.Robj")
gwa<-read.table("../data/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz",he=T)
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig")
m<-match(gwa$Snp,bim$V2)
gwa<-data.frame(gwa,bim[m,])
gwa$P<-gwa$MainP
gwa$SNP<-gwa$Snp
#mQTL<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr11.cistrans.txt.gz")
mQTL<-read.table(paste0("../data/",cpgid,".txt"),sep="\t",he=F)
h<-read.table("../data/header.txt",sep="\t",he=T)
names(mQTL)<-names(h)

regional_plot(diab_coloc=diab_coloc,snpid=snpid,cpgid=cpgid,trait=trait,gwa=gwa,mQTL=mQTL)
###

snpid="rs174528"
cpgid="cg16213375"
trait="Fasting_glucose"

load("../data/diab_coloc.Robj")
gwa<-read.table("../data/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz",he=T)
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig")
m<-match(gwa$Snp,bim$V2)
gwa<-data.frame(gwa,bim[m,])
gwa$P<-gwa$MainP
gwa$SNP<-gwa$Snp
#mQTL<-read.table("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.chr11.cistrans.txt.gz")
mQTL<-read.table(paste0("../data/",cpgid,".txt"),sep="\t",he=F)
h<-read.table("../data/header.txt",sep="\t",he=T)
names(mQTL)<-names(h)

regional_plot(diab_coloc=diab_coloc,snpid=snpid,cpgid=cpgid,trait=trait,gwa=gwa,mQTL=mQTL)
###
snpid="rs481887"
cpgid="cg23106115"
trait="Type 2 diabetes"


load("../data/diab_coloc.Robj")
gwa<-read.table("../data/diagram.mega-meta.txt",he=T)
gwa$chr<-gwa$CHR
bim<-read.table("/panfs/panasas01/shared-godmc/1kg_reference_ph3/eur.bim.orig")
m<-match(gwa$SNP,bim$V2)
gwa<-data.frame(gwa,bim[m,])


mQTL<-read.table(paste0("../data/",cpgid,".txt"),sep="\t",he=F)
h<-read.table("../data/header.txt",sep="\t",he=T)
names(mQTL)<-names(h)
regional_plot(diab_coloc=diab_coloc,snpid=snpid,cpgid=cpgid,trait=trait,gwa=gwa,mQTL=mQTL)

##




