#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)

r<-read.table(paste0(path,"selection_no_mhc_lct.txt"),he=T)
r$Type<-gsub("mqtls","All",r$Type)
cis<-read.table(paste0(path,"selection_cis_no_mhc_lct.txt"),he=T)
cis$Type<-gsub("cismqtls","cis only",cis$Type)
trans<-read.table(paste0(path,"selection_trans_no_mhc_lct.txt"),he=T)
trans$Type<-gsub("transmqtls","trans only",trans$Type)
amb<-read.table(paste0(path,"selection_ambivalent_no_mhc_lct.txt"),he=T)
amb$Type<-gsub("ambivalentmqtls","cis+trans",amb$Type)
cisall<-read.table(paste0(path,"selection_cis_no_lct_mhc_all.txt"),he=T)
cisall$Type<-gsub("cismqtls","cis any",cisall$Type)
transall<-read.table(paste0(path,"selection_trans_no_lct_mhc_all.txt"),he=T)
transall$Type<-gsub("transmqtls","trans any",transall$Type)

df<-rbind(r,cis,trans,amb,cisall,transall)
df$logOddsRatio<-log(df$OR)
df2<-df[df$PThresh=="1e-14",]
df2$Type<-as.character(df2$Type)

#w<-which(df2$Type=="mqtls")
#df2$Type[w]<-"All"
#df2$Type<-gsub("mqtls","",df2$Type)
#df2$Type<-gsub("cis","cis only",df2$Type)
#df2$Type<-gsub("trans","trans only",df2$Type)
#df2$Type<-gsub("ambivalent","cis+trans",df2$Type)

df2$Type<-as.factor(df2$Type)

pval_lim<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/mqtl_sds/garfield.Meff.mqtl_sds.out")
pval_lim<-pval_lim$V2[2]

p1<-ggplot(df2,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=Annotation))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
  ggsave(p1,file="./images/selection_no_mhc_lct.pdf",height=6,width=6)

pl3<-ggplot(df2,aes(x=Annotation,y=-log10(Pvalue),size=logOddsRatio,fill=Type))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
  ggsave(pl3,file="./images/selection2_no_mhc_lct.pdf",height=6,width=8)

pl3<-ggplot(df2,aes(x=Annotation,y=OR,size=-log10(Pvalue),fill=Type))+
  geom_hline(yintercept=1,col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,6),limits=c(0.2,6)) +
  ylab("Odds ratio (log scale)") +
  theme(legend.text=element_text(size=12))
  ggsave(pl3,file="./images/selection2_OR_no_mhc_lct.pdf",height=6,width=8)

  df$logOddsRatio<-log(df$OR)
  
#df3<-df2[which(df2$Pvalue<pval_lim),c("Annotation","OR","Pvalue","Beta","SE","CI95_lower","CI95_upper","NAnnotThesh","NAnnot","NThresh","N","Type")]
df3<-df2[,c("Annotation","OR","Pvalue","Beta","SE","CI95_lower","CI95_upper","NAnnotThesh","NAnnot","NThresh","N","Type")]
names(df3)<-gsub("NAnnotThesh","NAnnotThresh",names(df3))
names(df3)<-gsub("Annotation","Selection_metric",names(df3))
names(df3)<-gsub("Type","Annotation",names(df3))
df3$Selection_metric<-gsub("xpehhyri","XP-EHH(YRI)",df3$Selection_metric)
df3$Selection_metric<-gsub("xpehhchb","XP-EHH(CHB)",df3$Selection_metric)
df3$Selection_metric<-gsub("ihs","iHS",df3$Selection_metric)
df3$Selection_metric<-gsub("sds","SDS",df3$Selection_metric)
df3$Selection_metric<-gsub("fst","Fst",df3$Selection_metric)

w<-which(df3$Annotation=="All")
write.table(df3,"selection_results_nomhc_lct.txt",sep="\t",quote=F,row.names=F,col.names=T)


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

#hla<-which(f.all$snpchr=="chr6"&f.all$min>29570005&f.all$max<33377657)
#lct<-which(f.all$snpchr=="chr2"&f.all$min>134608646&f.all$max<138608646)

hla<-which(f.all$snpchr=="chr6"&f.all$min>24570005&f.all$max<38377657)
lct<-which(f.all$snpchr=="chr2"&f.all$min>129608646&f.all$max<143608646)

w<-c(hla,lct)
f.all<-f.all[-w,]

w<-which(f.all$mqtl_clumped=="TRUE")
f.all2<-f.all[w,]

m<-max(-log10(df2$Pvalue))

f.all2_cis<-f.all2[f.all2$Category=="cisonly",]

f.all2_sds<-f.all2[which(f.all2$sds=="1"),]


p1<-ggplot(f.all2, aes(x=sds_score, y=abs(max_abs_Effect))) +
  stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE,n=200) + 
  geom_point(size=0.5) +
  facet_wrap(~Category,scale="free_x",ncol=1)+
  stat_density2d(geom="tile", aes(fill=..density..^0.25,alpha=ifelse(..density..^0.25<0.4,0,1)), contour=FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))+
  xlim(-10,10)
  ggsave(p1,file="./images/selection_smoothscatter.pdf")

#dat_sig <- subset(dat, !grepl("metabolites__", fn) & nsnp > 3 & binom4 < 0.05/nrow(dat))
#dat_nsig <- subset(dat, !grepl("metabolites__", fn) & nsnp > 3 & binom4 >= 0.05/nrow(dat))

#p1 <- ggplot(subset(dat_nsig, !grepl("metabolites__", fn) & nsnp > 3), aes(x=label, y=-log10(binom4))) +
#geom_point(aes(size=nsnp)) +
#geom_point(data=dat_sig, aes(colour=as.factor(clust), size=nsnp)) +


w<-which(!is.na(f.all2$sds_score))
table(f.all2$Category[w])
#cisonly cis+trans transonly 
#    91837     25592       378 

p2<-ggplot(f.all2, aes(sds_score, abs(max_abs_Effect))) +
geom_point(alpha=0.2,size=0.5) +

#geom_point(data=hla,aes(colour="HLA"),size=0.5) +
#geom_point(data=lct,aes(colour="LCT"),size=0.5) +
#geom_point(data=wdfy4,aes(colour="WDFY4"),size=0.5) +

# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
facet_wrap(~Category,scale="free_x",ncol=3)+
geom_smooth(method="lm", colour="red") +
#geom_smooth(method="lm", formula=sds_score ~ MAF) +
labs(x="sds_score", y="max Effect size (SD)") +
xlim(-10,10)
ggsave(p2,file="./images/selection_sds_no_lct_mhc.png",height=2,width=8)

w<-which(!is.na(f.all2$Fst_score))
table(f.all2$Category[w])

#cisonly cis+trans transonly 
#   140579     54131       696 
p1<-ggplot(df2,aes(x=Annotation,y=logOddsRatio,size=-log10(Pvalue),fill=Type))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))

#
p2<-ggplot(f.all2, aes(sds_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour=Category)) +
#geom_smooth(method="lm", formula=sds_score ~ MAF) +

labs(x="sds_score", y="max Effect size")

p3<-ggplot(f.all2, aes(Fst_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour=Category)) +
#geom_smooth(method="lm", formula=Fst_score ~ MAF) +

labs(x="Fst_score", y="max Effect size")

ihs<-which(f.all2$iHS_score<5)
p4<-ggplot(f.all2[ihs,], aes(iHS_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour=Category)) +
#geom_smooth(data=f.all2,method="lm", formula=iHS_score ~ MAF) +

labs(x="iHS_score", y="max Effect size")

p5<-ggplot(f.all2, aes(xpehhchb_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="xpehhchb_score", y="max Effect size")

p6<-ggplot(f.all2, aes(xpehhyri_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="xpehhyri_score", y="max Effect size")

p7<-ggplot(f.all2, aes(Fst_score, MAF,colour=Category)) +
geom_point(alpha=0.1) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="Fst_score", y="MAF")

#png("../images/selection.png",height=1000,width=1000)
#grid.arrange(p1,p2,p3, layout_matrix = rbind(c(1,1),c(2,3)))
#g<-arrangeGrob(p1,p2,p3,layout_matrix = rbind(c(1,1),c(2,3)))
#ggsave(file="../images/selection.png",height=1000,width=1000)
#dev.off()

#png("../images/selection.png",height=1000,width=1800,res=300)
#pdf("../images/selection.pdf",height=30,width=18)
pdf("./images/selection.pdf",height=30,width=18, units="in", res=300)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1:2))
print(p2, vp = vplayout(2, 1))
print(p3, vp = vplayout(2, 2))
print(p4, vp = vplayout(3, 1))
print(p5, vp = vplayout(3, 2))
print(p6, vp = vplayout(4, 1))
#print(p7, vp = vplayout(4, 1))
dev.off()

##
model <- lm(data = f.all2, abs(max_abs_Effect) ~ sds_score + MAF)
f.all2$model <- stats::predict(model, newdata=f.all2)
err <- stats::predict(model, newdata=f.all2, se = TRUE)
f.all2$ucl <- err$fit + 1.96 * err$se.fit
f.all2$lcl <- err$fit - 1.96 * err$se.fit

p2<-ggplot(f.all2, aes(sds_score, model,colour=Category)) +
geom_point(alpha=0.1) +
geom_smooth(data=f.all2, method="lm",aes(x=sds_score, y=model, colour = Category,  ymin=lcl, ymax=ucl), size = 1.5,  se = TRUE) +
labs(x="sds_score", y="max Effect size")
ggsave(p2,file="./images/selectionSDS_MAFadj.png")
###
p2<-ggplot(f.all2, aes(sds_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
geom_hex( bins=30 ) +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="sds_score", y="max Effect size")
ggsave(p2,file="./images/selectionSDS_hexbin.png")

p2<-ggplot(f.all2, aes(sds_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
geom_density2d() +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="sds_score", y="max Effect size")
ggsave(p2,file="./images/selectionSDS_2ddensity.png")


p2<-ggplot(f.all2, aes(sds_score, min_log10pval,colour=Category)) +
geom_point(alpha=0.4) +
stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_tile() +
#scale_y_continuous(expand = c(0,0)) +
#scale_fill_gradient2(limits=c(-0.04,0.04),midpoint=mean(f.all2$density)) +
labs(x="sds_score", y="max -log10 Pvalue")

p3<-ggplot(f.all2, aes(Fst_score, min_log10pval,colour=Category)) +
geom_point(alpha=0.4) +
stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_tile() +
#scale_y_continuous(expand = c(0,0)) +
#scale_fill_gradient2(limits=c(-0.04,0.04),midpoint=mean(f.all2$density)) +
labs(x="Fst_score", y="max -log10 Pvalue")

pdf("./images/selection_pval.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1:2))
print(p2, vp = vplayout(2, 1))
print(p3, vp = vplayout(2, 2))
dev.off()


w1<-which(f.all$sds_score>2)
w2<-which(f.all$sds_score<=-2)
mean(abs(f.all$max_abs_Effect[w1]),na.rm=T)
#[1] 0.3965969
mean(abs(f.all$max_abs_Effect[w2]),na.rm=T)
#0.3979109
mean(abs(f.all$max_abs_Effect[-w1]),na.rm=T)
#0.4391176
mean(abs(f.all$max_abs_Effect[-w2]),na.rm=T)
#0.4390677

w<-which(!is.na(f.all$max_abs_Effect))

min(f.all$sds_score[w],na.rm=T)
#[1] -5.152155
max(f.all$sds_score[w],na.rm=T)
#[1] 6.623955
min(f.all$sds_score,na.rm=T)
#[1] -5.152155
max(f.all$sds_score,na.rm=T)
#[1] 6.623955



summary(lm(data=f.all2,sds_score~MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.058980   0.006397   9.220   <2e-16 ***
#MAF         -0.185387   0.021973  -8.437   <2e-16 ***

summary(lm(data=f.all2,abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.5836084  0.0016503 353.639   <2e-16 ***
#sds_score   -0.0010908  0.0007513  -1.452    0.147    
#MAF         -0.3168627  0.0056680 -55.904   <2e-16 ***

summary(lm(data=f.all2,abs(max_abs_Effect)~sds_score+nproxies+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.509e-01  1.653e-03 333.367   <2e-16 ***
#sds_score   -6.705e-04  7.305e-04  -0.918    0.359    
#nproxies     1.351e-03  1.637e-05  82.551   <2e-16 ***
#MAF         -3.265e-01  5.512e-03 -59.230   <2e-16 ***

summary(lm(data=f.all2,abs(max_abs_Effect)~sds_score+nproxies+tssdist+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.509e-01  1.653e-03 333.342   <2e-16 ***
#sds_score   -6.719e-04  7.305e-04  -0.920    0.358    
#nproxies     1.351e-03  1.637e-05  82.552   <2e-16 ***
#tssdist      5.320e-09  5.760e-09   0.924    0.356    
#MAF         -3.265e-01  5.512e-03 -59.226   <2e-16 ***

summary(lm(data=f.all2,abs(max_abs_Effect)~sds_score+nproxies+tssdist+GC_freq+CpG_freq+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  3.299e-01  4.846e-03  68.077   <2e-16 ***
#sds_score    4.755e-05  7.022e-04   0.068   0.9460    
#nproxies     1.537e-03  1.587e-05  96.842   <2e-16 ***
#tssdist      1.157e-08  5.537e-09   2.090   0.0366 *  
#GC_freq      2.357e-01  1.200e-02  19.646   <2e-16 ***
#CpG_freq     8.463e+00  1.280e-01  66.098   <2e-16 ***
#MAF         -3.306e-01  5.299e-03 -62.395   <2e-16 ***


#cis
summary(lm(data=f.all2[w1,],abs(max_abs_Effect)~sds_score+nproxies+tssdist+GC_freq+CpG_freq+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  3.343e-01  5.117e-03  65.334   <2e-16 ***
#sds_score   -7.442e-04  7.873e-04  -0.945    0.345    
#nproxies     1.502e-03  1.991e-05  75.409   <2e-16 ***
#tssdist      9.437e-09  6.069e-09   1.555    0.120    
#GC_freq      2.162e-01  1.244e-02  17.379   <2e-16 ***
#CpG_freq     7.708e+00  1.390e-01  55.449   <2e-16 ***
#MAF         -3.493e-01  5.861e-03 -59.598   <2e-16 ***

summary(lm(data=f.all2,Fst_score~nproxies+tssdist+GC_freq+CpG_freq))

#Call:
#lm(formula = Fst_score ~ nproxies + tssdist + GC_freq + CpG_freq, 
#    data = f.all2)

#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.18090 -0.08665 -0.03523  0.05193  0.73049 

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  1.120e-01  1.735e-03  64.559   <2e-16 ***
#nproxies     8.548e-05  5.391e-06  15.857   <2e-16 ***
#tssdist     -4.102e-09  2.126e-09  -1.929   0.0537 .  
#GC_freq      1.096e-02  4.409e-03   2.487   0.0129 *  
#CpG_freq     7.742e-02  3.862e-02   2.005   0.0450 *  

w0<-1:nrow(f.all2)
w1<-which(f.all2$snp_cis=="TRUE") #155213
w2<-which(f.all2$snp_cis=="FALSE") #789
w3<-which(f.all2$snp_cis=="ambivalent") #58599


df.out<-data.frame()
df.out2<-data.frame()
df.out_nomaf<-data.frame()
df.out2_nomaf<-data.frame()
outcomes<-c("sds_score","Fst_score","iHS_score","xpehhchb_score","xpehhyri_score")
cis<-c("All","TRUE","ambivalent","FALSE")
for (i in 1:length(outcomes)){
w<-which(names(f.all2)%in%outcomes[i])
for (j in 1:length(cis)){
if(cis[j]!="All"){w2<-which(f.all2$snp_cis==cis[j])}
if(cis[j]=="All"){w2<-1:nrow(f.all2)}

lmfit<-summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~f.all2[w2,w]+nproxies+tssdist+GC_freq+CpG_freq))
df<-data.frame(outcome=names(f.all2)[w],cis=cis[j],lmfit$coefficients)
df.out<-rbind(df.out,df) #complete model
df.out2<-rbind(df.out2,df[2,])

lmfit2<-summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~f.all2[w2,w]+nproxies+tssdist+GC_freq+CpG_freq+MAF))
df_nomaf<-data.frame(outcome=names(f.all2)[w],cis=cis[j],lmfit2$coefficients)
df.out_nomaf<-rbind(df.out_nomaf,df_nomaf)
df.out2_nomaf<-rbind(df.out2_nomaf,df_nomaf[2,])

}
}

df.out_maf<-data.frame()
df.out2_maf<-data.frame()
for (j in 1:length(cis)){
if(cis[j]!="All"){w2<-which(f.all2$snp_cis==cis[j])}
if(cis[j]=="All"){w2<-1:nrow(f.all2)}

lmfit3<-summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~f.all2[w2,"MAF"]+nproxies+tssdist+GC_freq+CpG_freq))
df_maf<-data.frame(cis=cis[j],lmfit3$coefficients)
df.out_maf<-rbind(df.out_maf,df_maf)
df.out2_maf<-rbind(df.out2_maf,df_maf[2,])
}
df.out3_maf<-df.out2_maf[which(df.out2_maf$cis!="All"),] 

#
w<-which(abs(f.all2$sds_score)>4)
f.all2[w,"SNP"]

w<-which(is.na(f.all2$sds_score))
f.all2$sds_res<-NULL
f.all2$fst_res<-NULL
f.all2$ihs_res<-NULL
f.all2$xpehhchb_res<-NULL
f.all2$xpehhyri_res<-NULL
f.all2$sds_res_nomaf<-NULL
f.all2$fst_res_nomaf<-NULL
f.all2$ihs_res_nomaf<-NULL
f.all2$xpehhchb_res_nomaf<-NULL
f.all2$xpehhyri_res_nomaf<-NULL


test<-data.frame(f.all2$sds_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2[a,],sds_score~nproxies+tssdist+GC_freq+CpG_freq+MAF))
f.all2[a,c("sds_res")]<-y

test<-data.frame(f.all2$Fst_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,Fst_score~nproxies+tssdist+GC_freq+CpG_freq+MAF))
f.all2[a,c("fst_res")]<-y

test<-data.frame(f.all2$iHS_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,iHS_score~nproxies+tssdist+GC_freq+CpG_freq+MAF))
f.all2[a,c("ihs_res")]<-y

test<-data.frame(f.all2$xpehhchb_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,xpehhchb_score~nproxies+tssdist+GC_freq+CpG_freq+MAF))
f.all2[a,c("xpehhchb_res")]<-y

test<-data.frame(f.all2$xpehhyri_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,xpehhyri_score~nproxies+tssdist+GC_freq+CpG_freq+MAF))
f.all2[a,c("xpehhyri_res")]<-y

#unadjusted for MAF
df.out2$cis<-gsub("TRUE","cisonly",df.out2$cis)
df.out2$cis<-gsub("ambivalent","cis+trans",df.out2$cis)
df.out2$cis<-gsub("FALSE","transonly",df.out2$cis)
o<-order(paste(df.out2$outcome,df.out2$cis))
df.out2<-df.out2[o,]
df.out3<-df.out2[which(df.out2$cis!="All"),] 

labs<-signif(df.out3[c(7:9,1:3,4:6,10:15),c("Pr...t.."),3],digits=3)
len <- length(unique(df.out3$cis))
dat <- data.frame(x = c(3.5,3.5,3.5,0.8,0.8,0.8,4.2,4.2,4.2,3.5,3.5,3.5,4.5,4.5,4.5), y = rep(2.5, len), Category=c("cisonly","cis+trans","transonly"),labs=labs)

#adjusted for MAF
df.out2_nomaf$cis<-gsub("TRUE","cisonly",df.out2_nomaf$cis)
df.out2_nomaf$cis<-gsub("ambivalent","cis+trans",df.out2_nomaf$cis)
df.out2_nomaf$cis<-gsub("FALSE","transonly",df.out2_nomaf$cis)
o<-order(paste(df.out2_nomaf$outcome,df.out2_nomaf$cis))
df.out2_nomaf<-df.out2_nomaf[o,]
df.out3_nomaf<-df.out2_nomaf[which(df.out2_nomaf$cis!="All"),] 

labs_nomaf<-signif(df.out3_nomaf[c(7:9,1:3,4:6,10:15),c("Pr...t.."),3],digits=3)
len_nomaf <- length(unique(df.out3_nomaf$cis))
dat_nomaf <- data.frame(x = c(2.5,2.5,2.5,0.8,0.8,0.8,2.5,2.5,2.5,2.5,2.5,2.5,4,4,4), y = rep(2.5, len_nomaf), Category=c("cisonly","cis+trans","transonly"),labs=labs_nomaf)



s<-which(!is.na(f.all2$sds_score))
p2<-ggplot(f.all2[s,], aes(sds_res, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
facet_wrap(~ Category) +
geom_smooth(method="lm", aes(colour="black")) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs_nomaf[1:3]), group=NULL),data=dat_nomaf[1:3,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="SDS_score", y="max Effect size")

s<-which(!is.na(f.all2$Fst_score))
p3<-ggplot(f.all2[s,], aes(fst_res, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
#geom_smooth(method="lm", formula=Fst_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs_nomaf[4:6]), group=NULL),data=dat_nomaf[4:6,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+
scale_x_continuous(limits = c(-0.25, 1))+

labs(x="Fst_score", y="max Effect size")

ihs<-which(f.all2$iHS_score<5)
p4<-ggplot(f.all2[ihs,], aes(ihs_res, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
#geom_smooth(data=f.all2,method="lm", formula=iHS_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs_nomaf[7:9]), group=NULL),data=dat_nomaf[7:9,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="iHS_score", y="max Effect size")

s<-which(!is.na(f.all2$xpehhchb_score))
p5<-ggplot(f.all2[s,], aes(xpehhchb_res, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs_nomaf[10:12]), group=NULL),data=dat_nomaf[10:12,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="XP_EHH(CHB)_score", y="max Effect size")

s<-which(!is.na(f.all2$xpehhyri_score))
p6<-ggplot(f.all2[s,], aes(xpehhyri_res, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs_nomaf[13:15]), group=NULL),data=dat_nomaf[13:15,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+
labs(x="XP_EHH(YRI)_score", y="max Effect size")

s<-which(!is.na(f.all2$Fst_score))
p7<-ggplot(f.all2[s,], aes(fst_res, MAF,colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +

labs(x="Fst_score", y="MAF")

png("./images/selection_MAFadj.png",height=30,width=16, units="in", res=300)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#print(p1, vp = vplayout(1, 1:2))
print(p2, vp = vplayout(1, 1:2))
print(p3, vp = vplayout(2, 1:2))
print(p4, vp = vplayout(3, 1:2))
print(p5, vp = vplayout(4, 1:2))
print(p6, vp = vplayout(5, 1:2))
dev.off()

#not adjusted for MAF
f.all2$sds_res_nomaf<-NULL
f.all2$fst_res_nomaf<-NULL
f.all2$ihs_res_nomaf<-NULL
f.all2$xpehhchb_res_nomaf<-NULL
f.all2$xpehhyri_res_nomaf<-NULL

test<-data.frame(f.all2$sds_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,sds_score~nproxies+tssdist+GC_freq+CpG_freq))
f.all2[a,c("sds_res_nomaf")]<-y

test<-data.frame(f.all2$Fst_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,Fst_score~nproxies+tssdist+GC_freq+CpG_freq))
f.all2[a,c("fst_res_nomaf")]<-y

test<-data.frame(f.all2$iHS_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,iHS_score~nproxies+tssdist+GC_freq+CpG_freq))
f.all2[a,c("ihs_res_nomaf")]<-y

test<-data.frame(f.all2$xpehhchb_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,xpehhchb_score~nproxies+tssdist+GC_freq+CpG_freq))
f.all2[a,c("xpehhchb_res_nomaf")]<-y

test<-data.frame(f.all2$xpehhyri_score,f.all2$nproxies,f.all2$tssdist,f.all2$GC_freq,f.all2$CpG_freq,f.all2$MAF)
a<-which(rowSums(is.na(test)) == 0)
y<-residuals(lm(data=f.all2,xpehhyri_score~nproxies+tssdist+GC_freq+CpG_freq))
f.all2[a,c("xpehhyri_res_nomaf")]<-y

s<-which(!is.na(f.all2$sds_score))
p2<-ggplot(f.all2[s,], aes(sds_res_nomaf, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

#geom_smooth(method="lm", formula=sds_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[1:3,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="SDS_score", y="max Effect size")

s<-which(!is.na(f.all2$Fst_score))
p3<-ggplot(f.all2[s,], aes(fst_res_nomaf, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

#geom_smooth(method="lm", formula=Fst_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[4:6,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+
scale_x_continuous(limits = c(-0.25, 1))+

labs(x="Fst_score", y="max Effect size")

ihs<-which(f.all2$iHS_score<5)
p4<-ggplot(f.all2[ihs,], aes(ihs_res_nomaf, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

#geom_smooth(data=f.all2,method="lm", formula=iHS_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[7:9,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="iHS_score", y="max Effect size")

s<-which(!is.na(f.all2$xpehhchb_score))
p5<-ggplot(f.all2[s,], aes(xpehhchb_res_nomaf, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[10:12,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="XP_EHH(CHB)_score", y="max Effect size")

s<-which(!is.na(f.all2$xpehhyri_score))
p6<-ggplot(f.all2[s,], aes(xpehhyri_res_nomaf, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[13:15,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="XP_EHH(YRI)_score", y="max Effect size")

s<-which(!is.na(f.all2$Fst_score))
p7<-ggplot(f.all2[s,], aes(MAF, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", signif(df.out3_maf[,5],digits=3)), group=NULL),data=data.frame(x=0.35,y=2.5),size=8,color = "black") +
labs(x="MAF", y="max Effect size")
ggsave(p7, file="./images/test.pdf", width=10, height=10)

#png("./images/selection_noMAF.png",height=30,width=18, units="in", res=300)
#create layout, assign it to viewport, push viewport to plotting device
#pdf("./images/selection_noMAF.pdf",height=10,width=8)
png("./images/selection_noMAF.png",height=30,width=18, units="in", res=300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#print(p1, vp = vplayout(1, 1:2))
print(p2, vp = vplayout(1, 1:2))
print(p3, vp = vplayout(2, 1:2))
print(p4, vp = vplayout(3, 1:2))
print(p5, vp = vplayout(4, 1:2))
print(p6, vp = vplayout(5, 1:2))
print(p7, vp = vplayout(6, 1:2))
dev.off()

#outliers
mean(abs(f.all2$tssdist))
#[1] 45797.88
mean(abs(f.all2$nproxies))
#[1] 26.58884
mean(abs(f.all2$GC_freq))
#[1] 0.4560906
mean(abs(f.all2$CpG_freq))
#[1] 0.01520923
mean(f.all2$iHS_score,na.rm=T)
#[1] 0.7026665
mean(f.all2$xpehhyri_score,na.rm=T)
#-0.2011362

test<-(f.all2[which(f.all2$xpehhyri_res_nomaf<(-2.5)&abs(f.all2$max_abs_Effect)>1),])
#test<-test[,c("closestgene","iHS_score","iHS_pval","ihs","ihs_res_nomaf","tssdist","GC_freq","CpG_freq","nproxies")]
mean(abs(test$tssdist))
#[1] 19362.57
mean(test$nproxies)
#[1] 925.7297
mean(test$GC_freq)
#[1] 0.4757012
mean(test$CpG_freq)
#[1] 0.015904
mean(test$iHS_score,na.rm=T)
#[1] 0.491596
mean(test$xpehhyri_score,na.rm=T)
#-1.67363

s<-which(!is.na(f.all2$sds_score))
p2<-ggplot(f.all2[s,], aes(sds_res_nomaf, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

#geom_smooth(method="lm", formula=sds_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[1:3,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="SDS_score", y="max Effect size")

s<-which(!is.na(f.all2$Fst_score))
p3<-ggplot(f.all2[s,], aes(Fst_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

#geom_smooth(method="lm", formula=Fst_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[4:6,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+
scale_x_continuous(limits = c(-0.25, 1))+

labs(x="Fst_score", y="max Effect size")

ihs<-which(f.all2$iHS_score<5)
p4<-ggplot(f.all2[ihs,], aes(iHS_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

#geom_smooth(data=f.all2,method="lm", formula=iHS_score ~ MAF) +
theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[7:9,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="iHS_score", y="max Effect size")

s<-which(!is.na(f.all2$xpehhchb_score))
p5<-ggplot(f.all2[s,], aes(xpehhchb_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[10:12,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="XP_EHH(CHB)_score", y="max Effect size")

s<-which(!is.na(f.all2$xpehhyri_score))
p6<-ggplot(f.all2[s,], aes(xpehhyri_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", labs), group=NULL),data=dat[13:15,],size=8,color = "black") +
scale_y_continuous(limits = c(0, 2.5))+

labs(x="XP_EHH(YRI)_score", y="max Effect size")

s<-which(!is.na(f.all2$Fst_score))
p7<-ggplot(f.all2[s,], aes(MAF, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.1) +
facet_wrap(~ Category) +
# stat_density_2d(aes(fill = ..level..), geom="polygon") +
#geom_smooth(aes(colour=Category)) +
geom_smooth(method="lm", aes(colour="black")) +
geom_smooth(aes(colour="red")) +

theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.text.y = element_text( size = 20 ),
           axis.title = element_text( size = 20, face = "bold" ),
           legend.position="none",strip.text = element_text(size = 20)) +
geom_text(aes(x, y, label=paste("p =", signif(df.out3_maf[,5],digits=3)), group=NULL),data=data.frame(x=0.4,y=2.5),size=8,color = "black") +
labs(x="MAF", y="max Effect size")
ggsave(p7, file="./images/test.pdf", width=10, height=10)

#png("./images/selection_noMAF.png",height=30,width=18, units="in", res=300)
#create layout, assign it to viewport, push viewport to plotting device
#pdf("./images/selection_noMAF.pdf",height=10,width=8)
png("./images/selection_noMAFforpaper.png",height=30,width=18, units="in", res=300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#print(p1, vp = vplayout(1, 1:2))
print(p2, vp = vplayout(1, 1:2))
print(p3, vp = vplayout(2, 1:2))
print(p4, vp = vplayout(3, 1:2))
print(p5, vp = vplayout(4, 1:2))
print(p6, vp = vplayout(5, 1:2))
print(p7, vp = vplayout(6, 1:2))
dev.off()

