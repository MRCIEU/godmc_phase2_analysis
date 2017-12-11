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


pval_lim<-0.01
p1<-ggplot(df2,aes(x=Type,y=-log10(Pvalue),size=logOddsRatio,fill=Annotation))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
  ggsave(p1,file="selection.pdf",height=6,width=6)

pl3<-ggplot(df2,aes(x=Annotation,y=-log10(Pvalue),size=logOddsRatio,fill=Type))+
  geom_hline(yintercept=-log10(pval_lim),col="black",linetype="dashed")+
  geom_point(alpha=0.7,shape=21,stroke=1)+
  facet_wrap(~Category,scale="free_x",ncol=1)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="bottom")+
  scale_size(range=c(1,4))+
  scale_color_manual(values=c("TRUE"="black","FALSE"="#EEEEEE"))
  ggsave(pl3,file="selection2.pdf",height=6,width=6)

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


f.all2<-f.all[-w,]

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


#png("../images/selection.png",height=1000,width=1000)
#grid.arrange(p1,p2,p3, layout_matrix = rbind(c(1,1),c(2,3)))
#g<-arrangeGrob(p1,p2,p3,layout_matrix = rbind(c(1,1),c(2,3)))
#ggsave(file="../images/selection.png",height=1000,width=1000)
#dev.off()

#png("../images/selection.png",height=1000,width=1800,res=300)
#pdf("../images/selection.pdf",height=30,width=18)
png("../images/selection.png",height=30,width=18, units="in", res=300)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1:2))
#print(p2, vp = vplayout(2, 1))
#print(p3, vp = vplayout(2, 2))
#print(p4, vp = vplayout(3, 1))
print(p5, vp = vplayout(3, 2))
#print(p6, vp = vplayout(4, 1))
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
ggsave(p2,file="../images/selectionSDS_MAFadj.png")
###
p2<-ggplot(f.all2, aes(sds_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
geom_hex( bins=30 ) +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="sds_score", y="max Effect size")
ggsave(p2,file="../images/selectionSDS_hexbin.png")

p2<-ggplot(f.all2, aes(sds_score, abs(max_abs_Effect),colour=Category)) +
geom_point(alpha=0.2) +
geom_density2d() +
geom_smooth(method="lm", aes(colour=Category)) +
labs(x="sds_score", y="max Effect size")
ggsave(p2,file="../images/selectionSDS_2ddensity.png")




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

pdf("../images/selection_pval.pdf",height=10,width=18)
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
#[1] 0.2191319
mean(abs(f.all$max_abs_Effect[w2]),na.rm=T)
#0.2126606
mean(abs(f.all$max_abs_Effect[-w1]),na.rm=T)
#0.2573107
mean(abs(f.all$max_abs_Effect[-w2]),na.rm=T)
#0.2572907

w<-which(!is.na(f.all$max_abs_Effect))

min(f.all$sds_score[w],na.rm=T)
#[1] -7.433762
max(f.all$sds_score[w],na.rm=T)
#[1] 8.568006
min(f.all$sds_score,na.rm=T)
#[1] -9.241373
max(f.all$sds_score,na.rm=T)
#[1] 10.00023


#minus HLA
#chr6:29570005-33377657
f.all$sds_scoresq<-f.all$sds_score^2
hla<-which(f.all$snpchr=="chr6"&f.all$min>29570005&f.all$max<33377657)
lct<-which(f.all$snpchr=="chr2"&f.all$min>134608646&f.all$max<138608646)
excl<-c(hla,lct)
f.all3<-f.all[-excl,]

w1<-which(f.all$snp_cis=="TRUE") #216244
w2<-which(f.all$snp_cis=="FALSE") #13516
w3<-which(f.all$snp_cis=="ambivalent") #2842
w4<-which(f.all$snp_cis!="TRUE") #16358

w1_3<-which(f.all3$snp_cis=="TRUE") #216244
w2_3<-which(f.all3$snp_cis=="FALSE") #13516
w3_3<-which(f.all3$snp_cis=="ambivalent") #2842
w4_3<-which(f.all3$snp_cis!="TRUE") #16358

summary(lm(data=f.all,sds_score~MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.0374604  0.0009767   38.35   <2e-16 ***
#MAF         -0.1549401  0.0035953  -43.09   <2e-16 ***

summary(lm(data=f.all,abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3053206  0.0009761 312.791  < 2e-16 ***
#sds_score   -0.0014383  0.0004248  -3.386  0.00071 ***
#MAF         -0.3287298  0.0034009 -96.661  < 2e-16 ***

summary(lm(formula = abs(max_abs_Effect) ~ sds_score, data = f.all))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.2211922  0.0004584 482.525   <2e-16 ***
#sds_score   0.0002190  0.0004403   0.497    0.619    

summary(lm(data=f.all,abs(max_abs_Effect)~sds_scoresq+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3067697  0.0010089 304.055  < 2e-16 ***
#sds_scoresq -0.0012761  0.0002074  -6.152 7.67e-10 ***
#MAF         -0.3290916  0.0034004 -96.781  < 2e-16 ***


summary(lm(data=f.all3,abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3051565  0.0009788 311.764  < 2e-16 ***
#sds_score   -0.0015236  0.0004311  -3.535 0.000409 ***
#MAF         -0.3283422  0.0034086 -96.328  < 2e-16 ***

summary(lm(data=f.all3,abs(max_abs_Effect)~sds_scoresq+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3070114  0.0010145 302.621  < 2e-16 ***
#sds_scoresq -0.0016529  0.0002234  -7.397  1.4e-13 ***
#MAF         -0.3288663  0.0034080 -96.497  < 2e-16 ***


#cis
summary(lm(data=f.all[w1,],abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3058840  0.0010144 301.553   <2e-16 ***
#sds_score   -0.0009688  0.0004534  -2.137   0.0326 *  
#MAF         -0.3299891  0.0035412 -93.185   <2e-16 ***

#trans
summary(lm(data=f.all[w2,],abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.279010   0.003309  84.323  < 2e-16 ***
#sds_score   -0.006256   0.001104  -5.667 1.51e-08 ***
#MAF         -0.304603   0.011178 -27.249  < 2e-16 ***

summary(lm(data=f.all3[w2,],abs(max_abs_Effect)~sds_score))

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.198030   0.001538 128.738   <2e-16 ***
#sds_score   -0.001959   0.001155  -1.696     0.09 .  

summary(lm(data=f.all[w2,],abs(max_abs_Effect)~abs(sds_score)+MAF))
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.296527   0.003727   79.57   <2e-16 ***
#abs(sds_score) -0.017118   0.001527  -11.21   <2e-16 ***
#MAF            -0.313547   0.011099  -28.25   <2e-16 ***

summary(lm(data=f.all3[w2_3,],abs(max_abs_Effect)~sds_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.279899   0.003323  84.241  < 2e-16 ***
#sds_score   -0.006105   0.001121  -5.446 5.34e-08 ***
#MAF         -0.306140   0.011219 -27.288  < 2e-16 ***



#ambivalent
summary(lm(data=f.all[w3,],abs(max_abs_Effect)~sds_score+MAF))

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.415255   0.011873  34.974   <2e-16 ***
#sds_score   -0.003960   0.003081  -1.285    0.199    
#MAF         -0.386760   0.039364  -9.825   <2e-16 ***

summary(lm(data=f.all[w3,],abs(max_abs_Effect)~sds_score))
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.309747   0.005287  58.588   <2e-16 ***
#sds_score   0.002597   0.003140   0.827    0.408    


#ambivalent+trans
summary(lm(data=f.all[w4,],abs(max_abs_Effect)~sds_score+MAF))

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.297237   0.003474  85.561  < 2e-16 ***
#sds_score   -0.004940   0.001109  -4.455  8.5e-06 ***
#MAF         -0.311027   0.011704 -26.575  < 2e-16 ***

summary(lm(data=f.all[w4,],abs(max_abs_Effect)~sds_score))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.2142367  0.0015923 134.545   <2e-16 ***
#sds_score   -0.0003745  0.0011467  -0.327    0.744    

##FST

summary(lm(data=f.all,abs(max_abs_Effect)~Fst_score+MAF))

#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3914197  0.0008361  468.15   <2e-16 ***
#Fst_score   -0.0535729  0.0037570  -14.26   <2e-16 ***
#MAF         -0.5619278  0.0030060 -186.94   <2e-16 ***

summary(lm(data=f.all3,abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3905328  0.0008479  460.56   <2e-16 ***
#Fst_score   -0.0546688  0.0037754  -14.48   <2e-16 ***
#MAF         -0.5569243  0.0030353 -183.48   <2e-16 ***

#cis
summary(lm(data=f.all[w1,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3931246  0.0008739  449.86   <2e-16 ***
#Fst_score   -0.0558333  0.0039460  -14.15   <2e-16 ***
#MAF         -0.5674103  0.0031523 -180.00   <2e-16 ***

summary(lm(data=f.all3[w1_3,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.3916384  0.0008844  442.85   <2e-16 ***
#Fst_score   -0.0553041  0.0039612  -13.96   <2e-16 ***
#MAF         -0.5620373  0.0031779 -176.86   <2e-16 ***


#trans
summary(lm(data=f.all[w2,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.345260   0.002780 124.216  < 2e-16 ***
#Fst_score   -0.035088   0.011932  -2.941  0.00328 ** 
#MAF         -0.481530   0.009662 -49.840  < 2e-16 ***

summary(lm(data=f.all[w2,],abs(max_abs_Effect)~Fst_score))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.244445   0.002087  117.14   <2e-16 ***
#Fst_score   -0.146972   0.012828  -11.46   <2e-16 ***

#ambivalent
summary(lm(data=f.all[w3,],abs(max_abs_Effect)~Fst_score+MAF))
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.490264   0.008018  61.144   <2e-16 ***
#Fst_score   -0.015418   0.031001  -0.497    0.619    
#MAF         -0.574796   0.026994 -21.293   <2e-16 ***

summary(lm(data=f.all[w3,],abs(max_abs_Effect)~Fst_score))
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
ggsave(p2, file="sdspval.pdf", width=10, height=10)

#




