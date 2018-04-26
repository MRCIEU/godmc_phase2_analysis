#/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/07_enrichments

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output/"
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)

r<-read.table(paste0(path,"selection.txt"),he=T)
cis<-read.table(paste0(path,"selection_cis.txt"),he=T)
trans<-read.table(paste0(path,"selection_trans.txt"),he=T)
amb<-read.table(paste0(path,"selection_ambivalent.txt"),he=T)

df<-rbind(r,cis,trans,amb)
df$logOddsRatio<-log(df$OR)
df2<-df[df$PThresh=="1e-14",]
df2$Type<-as.character(df2$Type)
w<-which(df2$Type=="mqtls")
df2$Type[w]<-"All"
df2$Type<-gsub("mqtls","",df2$Type)
df2$Type<-gsub("cis","cis only",df2$Type)
df2$Type<-gsub("trans","trans only",df2$Type)
df2$Type<-gsub("ambivalent","cis+trans",df2$Type)

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
  ggsave(p1,file="./images/selection.pdf",height=6,width=6)

pl3<-ggplot(df2,aes(x=Annotation,y=-log10(Pvalue),size=logOddsRatio,fill=Type))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
  ggsave(pl3,file="./images/selection2.pdf",height=6,width=8)

pl3<-ggplot(df2,aes(x=Annotation,y=OR,size=-log10(Pvalue),fill=Type))+
  geom_hline(yintercept=1,col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1,2,3,4,5,6),limits=c(0.2,6)) +
  ylab("Odds ratio (log scale)") +
  theme(legend.text=element_text(size=12))
  ggsave(pl3,file="./images/selection2_OR.pdf",height=6,width=8)

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

m<-max(-log10(df2$Pvalue))

p1<-ggplot(df2,aes(x=Annotation,y=-log10(Pvalue),size=logOddsRatio,fill=Type))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  ylim(0,(m+3))+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=15),legend.position="bottom")+
  scale_size(range=c(1,8))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))

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
png("./images/selection.png",height=30,width=18, units="in", res=300)
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
print(p7, vp = vplayout(4, 1))
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
#[1] 0.4048265
mean(abs(f.all$max_abs_Effect[w2]),na.rm=T)
#0.4054604
mean(abs(f.all$max_abs_Effect[-w1]),na.rm=T)
#0.4440296
mean(abs(f.all$max_abs_Effect[-w2]),na.rm=T)
#0.4439892

w<-which(!is.na(f.all$max_abs_Effect))

min(f.all$sds_score[w],na.rm=T)
#[1] -9.241373
max(f.all$sds_score[w],na.rm=T)
#[1] 10.00023
min(f.all$sds_score,na.rm=T)
#[1] -9.241373
max(f.all$sds_score,na.rm=T)
#[1] 10.00023


#minus HLA
#chr6:29570005-33377657
f.all2$sds_scoresq<-f.all2$sds_score^2
hla<-which(f.all2$snpchr=="chr6"&f.all2$min>29570005&f.all2$max<33377657)
lct<-which(f.all2$snpchr=="chr2"&f.all2$min>134608646&f.all2$max<138608646)
excl<-c(hla,lct)
f.all3<-f.all2[-excl,]
f.all4<-f.all2[-hla,]

w1<-which(f.all2$snp_cis=="TRUE") #5138085
w2<-which(f.all2$snp_cis=="FALSE") #95610
w3<-which(f.all2$snp_cis=="ambivalent") #971913
w4<-which(f.all2$snp_cis!="TRUE") #1067523

w1_3<-which(f.all3$snp_cis=="TRUE") #216244
w2_3<-which(f.all3$snp_cis=="FALSE") #13516
w3_3<-which(f.all3$snp_cis=="ambivalent") #2842
w4_3<-which(f.all3$snp_cis!="TRUE") #16358

w1_4<-which(f.all4$snp_cis=="TRUE") #216244
w2_4<-which(f.all4$snp_cis=="FALSE") #13516
w3_4<-which(f.all4$snp_cis=="ambivalent") #2842
w4_4<-which(f.all4$snp_cis!="TRUE") #16358


summary(lm(data=f.all2,sds_score~MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.093089   0.006607   14.09   <2e-16 ***
#MAF         -0.290297   0.022727  -12.77   <2e-16 ***

summary(lm(data=f.all2,abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.5910515  0.0016448 359.352  < 2e-16 ***
#sds_score    0.0026694  0.0007181   3.718 0.000201 ***
#MAF         -0.3318634  0.0056571 -58.663  < 2e-16 ***

summary(lm(formula = abs(max_abs_Effect) ~ sds_score, data = f.all2))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.5050888  0.0007576 666.686  < 2e-16 ***
#sds_score   0.0042215  0.0007278   5.801 6.63e-09 ***

summary(lm(data=f.all2,abs(max_abs_Effect)~sds_scoresq+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.5786976  0.0016957  341.28   <2e-16 ***
#sds_scoresq  0.0101624  0.0003539   28.72   <2e-16 ***
#MAF         -0.3264990  0.0056383  -57.91   <2e-16 ***

summary(lm(data=f.all3,abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.5890487  0.0016450 358.078  < 2e-16 ***
#sds_score    0.0024532  0.0007265   3.377 0.000733 ***
#MAF         -0.3292750  0.0056566 -58.211  < 2e-16 ***

summary(lm(data=f.all3,abs(max_abs_Effect)~sds_scoresq+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.5781901  0.0017017  339.77   <2e-16 ***
#sds_scoresq  0.0091772  0.0003769   24.35   <2e-16 ***
#MAF         -0.3246689  0.0056435  -57.53   <2e-16 ***

#cis
summary(lm(data=f.all2[w1,],abs(max_abs_Effect)~sds_score))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.4752920  0.0008303  572.41   <2e-16 ***
#sds_score   -0.0004456  0.0008400   -0.53    0.596    

summary(lm(data=f.all3[w1_3,],abs(max_abs_Effect)~sds_score))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.5015556  0.0009525 526.579  < 2e-16 ***
#sds_score   0.0046584  0.0009208   5.059 4.22e-07 ***

#trans
summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.327676   0.017799  18.410  < 2e-16 ***
#sds_score    0.003451   0.009475   0.364 0.715894    
#MAF         -0.219694   0.061514  -3.571 0.000401 ***

summary(lm(data=f.all3[w2_3,],abs(max_abs_Effect)~sds_score))

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.50365    0.01346  37.406   <2e-16 ***
#sds_score   -0.02500    0.01213  -2.062     0.04 *  

summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~abs(sds_score)+MAF))
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.35735    0.02034  17.573  < 2e-16 ***
#abs(sds_score) -0.04363    0.01496  -2.917 0.003747 ** 
#MAF            -0.21756    0.06084  -3.576 0.000394 ***

summary(lm(data=f.all3[w2_3,],abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.59248    0.02828  20.947  < 2e-16 ***
#sds_score   -0.02671    0.01193  -2.238 0.025848 *  
#MAF         -0.34705    0.09766  -3.554 0.000435 ***

#ambivalent
summary(lm(data=f.all2[w3,],abs(max_abs_Effect)~sds_score+MAF))

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.722237   0.003739 193.177  < 2e-16 ***
#sds_score    0.008596   0.001328   6.473 9.78e-11 ***
#MAF         -0.406790   0.012308 -33.052  < 2e-16 ***

summary(lm(data=f.all3[w3_3,],abs(max_abs_Effect)~sds_score))
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.504915   0.002243 225.126   <2e-16 ***
#sds_score   0.005207   0.002185   2.383   0.0172 *  


#ambivalent+trans
summary(lm(data=f.all2[w4,],abs(max_abs_Effect)~sds_score+MAF))

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.715127   0.003736 191.435  < 2e-16 ***
#sds_score    0.008739   0.001333   6.555 5.68e-11 ***
#MAF         -0.398454   0.012305 -32.381  < 2e-16 ***

summary(lm(data=f.all2[w4,],abs(max_abs_Effect)~sds_score))
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.605781   0.001628 372.164   <2e-16 ***
#sds_score   0.012196   0.001354   9.007   <2e-16 ***

##FST

summary(lm(data=f.all2,abs(max_abs_Effect)~Fst_score+MAF))

#Coefficients:
#             Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)  0.712370   0.001279  556.758  < 2e-16 ***
#Fst_score   -0.024809   0.005718   -4.339 1.43e-05 ***
#MAF         -0.576544   0.004515 -127.688  < 2e-16 ***

summary(lm(data=f.all3,abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.693635   0.001282  541.03   <2e-16 ***
#Fst_score    0.012660   0.005678    2.23   0.0258 *  
#MAF         -0.554447   0.004509 -122.96   <2e-16 ***

#cis
summary(lm(data=f.all2[w1,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)  0.658781   0.001416  465.085   <2e-16 ***
#Fst_score    0.012094   0.006534    1.851   0.0642 .  
#MAF         -0.613491   0.005137 -119.435   <2e-16 ***

summary(lm(data=f.all3[w1_3,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.691121   0.001616 427.613   <2e-16 ***
#Fst_score    0.020227   0.007155   2.827   0.0047 ** 
#MAF         -0.555094   0.005692 -97.517   <2e-16 ***


#trans
summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.345260   0.002780 124.216  < 2e-16 ***
#Fst_score   -0.035088   0.011932  -2.941  0.00328 ** 
#MAF         -0.481530   0.009662 -49.840  < 2e-16 ***

summary(lm(data=f.all2[w2,],abs(max_abs_Effect)~Fst_score))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.244445   0.002087  117.14   <2e-16 ***
#Fst_score   -0.146972   0.012828  -11.46   <2e-16 ***

#ambivalent
summary(lm(data=f.all2[w3,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.490264   0.008018  61.144   <2e-16 ***
#Fst_score   -0.015418   0.031001  -0.497    0.619    
#MAF         -0.574796   0.026994 -21.293   <2e-16 ***

summary(lm(data=f.all2[w3,],abs(max_abs_Effect)~Fst_score))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.361467   0.005687  63.559  < 2e-16 ***
#Fst_score   -0.132892   0.032960  -4.032 5.68e-05 ***

summary(lm(data=f.all,abs(max_abs_Effect)~xpehhchb_score+MAF))
#Coefficients:
#                 Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)     0.3879425  0.0008065  481.026  < 2e-16 ***
#xpehhchb_score -0.0022740  0.0004535   -5.015 5.31e-07 ***
#MAF            -0.5750004  0.0030703 -187.280  < 2e-16 ***

summary(lm(data=f.all,abs(max_abs_Effect)~xpehhyri_score+MAF))

#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.3883473  0.0008056  482.06   <2e-16 ***
#xpehhyri_score -0.0062458  0.0004654  -13.42   <2e-16 ***
#MAF            -0.5808613  0.0030914 -187.90   <2e-16 ***

summary(lm(data=f.all,abs(max_abs_Effect)~iHS_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3886247  0.0009277  418.92   <2e-16 ***
#iHS_score   -0.0137831  0.0007132  -19.32   <2e-16 ***
#MAF         -0.5439665  0.0030267 -179.72   <2e-16 ***

ihs<-which(f.all$iHS_score<5)
summary(lm(data=f.all[ihs,],abs(max_abs_Effect)~iHS_score+MAF))

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3887154  0.0009284  418.70   <2e-16 ***
#iHS_score   -0.0139444  0.0007161  -19.47   <2e-16 ***
#MAF         -0.5438825  0.0030269 -179.69   <2e-16 ***


df$min_pval2<-df$min_pval
w<-which(df$min_pval==0)
df$min_pval2[w]<-min(df$min_pval[-w],na.rm=T)

p2<-ggplot(df, aes(sds_score, -log10(df$min_pval2))) +
geom_point() +
theme_minimal()
ggsave(p2, file="./images/sdspval.pdf", width=10, height=10)

#
w<-which(abs(f.all2$sds_score)>4)
f.all2[w,"SNP"]


