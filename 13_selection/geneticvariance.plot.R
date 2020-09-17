library(ggplot2)

load("./height_sds/extremeheight_plot.Robj")
spl<-strsplit(as.character(h_all$name),":")
spl<-do.call("rbind",spl)
w1<-which(spl[,1]=="chr6"&spl[,2]>24570005&spl[,2]<38377657)
w2<-which(spl[,1]=="chr2"&spl[,2]>129608646&spl[,2]<143608646)
w<-(unique(w1,w2))
if(length(w)>0){
  h_all<-h_all[-w,]  
}
df1<-data.frame(SNP=h_all$name,trait="Extreme Height",es=h_all$es,mqtl=h_all$height_mqtl)
summary(lm(df1$es~df1$mqtl))

#Coefficients:
#                                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                           0.004884   0.001086   4.496 1.75e-05 ***
#df1$mqtlextreme height mqtl (n=48)    0.010198   0.001622   6.288 6.97e-09 ***
#df1$mqtlextreme height mqtl sds (n=4) 0.008280   0.004311   1.921   0.0574 .  

#Coefficients:
#                                        Estimate Std. Error t value Pr(>|t|)
#(Intercept)                            0.0133058  0.0007134  18.652   <2e-16
#df1$mqtlextreme height mqtl (n=48)     0.0019260  0.0011117   1.732   0.0854
#df1$mqtlextreme height mqtl sds (n=4) -0.0001426  0.0032690  -0.044   0.9653



load("./height_sds/height_Yenko_plot.Robj")
w1<-which(h_all$CHR=="6"&h_all$POS>24570005&h_all$POS<38377657)
w2<-which(h_all$CHR=="2"&h_all$POS>129608646&h_all$POS<143608646)
w<-(unique(w1,w2))
if(length(w)>0){
h_all<-h_all[-w,]  
}
df2<-data.frame(SNP=h_all$ID,trait="Height",es=h_all$es,mqtl=h_all$height_mqtl)
summary(lm(df2$es~df2$mqtl))
#Coefficients:
#                                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                     1.065e-04  2.604e-06  40.911  < 2e-16 ***
#df2$mqtlheight mQTL(n=3948)    -1.489e-05  3.511e-06  -4.240 2.27e-05 ***
#df2$mqtlheight mQTL SDS (n=44)  5.616e-05  2.354e-05   2.386   0.0171 *  

#Coefficients:
#                                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                     1.065e-04  2.549e-06  41.796  < 2e-16 ***
#df2$mqtlheight mQTL(n=3948)    -1.515e-05  3.282e-06  -4.617 3.96e-06 ***
#df2$mqtlheight mQTL SDS (n=44)  5.616e-05  2.304e-05   2.437   0.0148 *  

w<-which(df2$mqtl=="height SNPs (n=3290)")
summary(lm(df2$es[-w]~df2$mqtl[-w]))

w<-which(df2$mqtl=="height mQTL SDS (n=44)")
summary(lm(df2$es[-w]~df2$mqtl[-w]))



load("./cvd_sds/cvd_plot.Robj")
w1<-which(h_all$chr=="6"&h_all$bp_hg19>24570005&h_all$bp_hg19<38377657)
w2<-which(h_all$chr=="2"&h_all$bp_hg19>129608646&h_all$bp_hg19<143608646)
w<-(unique(w1,w2))
if(length(w)>0){
h_all<-h_all[-w,]  
}
df3<-data.frame(SNP=h_all$SNP,trait="CVD",es=h_all$es,mqtl=h_all$cvd_mqtl)
summary(lm(df3$es~df3$mqtl))
#Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                0.0015891  0.0001574  10.098  < 2e-16 ***
#df3$mqtlcvd mqtl (n=15)    0.0002737  0.0005916   0.463    0.644    
#df3$mqtlcvd mqtl sds (n=2) 0.0088471  0.0015697   5.636 5.53e-08 ***

#Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                0.0015891  0.0001694   9.382  < 2e-16 ***
#df3$mqtlcvd mqtl (n=15)    0.0009666  0.0003626   2.666  0.00817 ** 
#df3$mqtlcvd mqtl sds (n=2) 0.0088471  0.0016896   5.236 3.47e-07 ***

#chr9:22098619:SNP   CVD 1.869237e-02 cvd mqtl sds (n=2) #rs2891168
#chr12:112059557:SNP   CVD 2.179914e-03 cvd mqtl sds (n=2) #rs11065979


load("./cd_sds/cd_plot.Robj")
w1<-which(h_all$CHR=="6"&h_all$BP>24570005&h_all$BP<38377657)
w2<-which(h_all$CHR=="2"&h_all$BP>129608646&h_all$BP<143608646)
w<-(unique(w1,w2))
if(length(w)>0){
h_all<-h_all[-w,]  
}
df4<-data.frame(SNP=h_all$SNP,trait="Crohn's Disease",es=h_all$es,mqtl=h_all$cd_mqtl)

load("./scz_sds/scz_plot.Robj")
w1<-which(h_all$CHR=="6"&h_all$BP>24570005&h_all$BP<38377657)
w2<-which(h_all$CHR=="2"&h_all$BP>129608646&h_all$BP<143608646)
w<-(unique(w1,w2))
if(length(w)>0){
h_all<-h_all[-w,]  
}
df5<-data.frame(SNP=h_all$ID,trait="Schizophrenia",es=h_all$es,mqtl=h_all$scz_mqtl)

df<-rbind(df1,df2,df3,df4,df5)

df$mqtl<-gsub("mqtl","mQTL",df$mqtl)
df$mqtl<-gsub("sds","SDS",df$mqtl)
df$mqtl<-gsub("scz","Schizophrenia",df$mqtl)
df$mqtl<-gsub("extreme height","Extreme Height",df$mqtl)
df$mqtl<-gsub("height","Height",df$mqtl)
df$mqtl<-gsub("cvd","CVD",df$mqtl)
df$mqtl<-gsub("cd","Crohn's Disease",df$mqtl)

df$mqtl2<-df$mqtl
df$mqtl2<-gsub("Extreme Height ","",df$mqtl2,fixed=T)
df$mqtl2<-gsub("Height ","",df$mqtl2,fixed=T)
df$mqtl2<-gsub("CVD ","",df$mqtl2,fixed=T)
df$mqtl2<-gsub("Crohn's Disease ","",df$mqtl2,fixed=T)
df$mqtl2<-gsub("Schizophrenia ","",df$mqtl2,fixed=T)
df$mqtl2<-gsub("SNPs","All SNPs",df$mqtl2,fixed=T)
df$mqtl2<-gsub(" mQTL SDS ","mQTL SDS SNPs",df$mqtl2,fixed=T)

df$mqtl2<-gsub("All SNPs (n=60)","All SNPs (n=80)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("All SNPs (n=159)","All SNPs (n=164)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("All SNPs (n=202)","All SNPs (n=197)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("All SNPs (n=3290)","All SNPs (n=3229)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("All SNPs (n=31)","All SNPs (n=47)",df$mqtl2,fixed=T)

df$mqtl2<-gsub("mQTL (n=15)","mQTL (n=55)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("mQTL (n=226)","mQTL (n=152)",df$mqtl2,fixed=T)
g<-grep("3948",df$mqtl2)
df$mqtl2[g]<-"mQTL (n=4910)"
df$mqtl2<-gsub("mQTL (n=48)","mQTL (n=56)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("mQTL (n=60)","mQTL (n=75)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("mQTL (n=15)","mQTL (n=55)",df$mqtl2,fixed=T)

df$mqtl2<-gsub("mQTL SDS (n=44)","mQTL SDS (n=40)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("mQTL SDS (n=9)","mQTL SDS (n=4)",df$mqtl2,fixed=T)
df$mqtl2<-gsub("mQTL SDS (n=3)","mQTL SDS (n=2)",df$mqtl2,fixed=T)

table(df$mqtl2)
#df$mqtl2<-factor(df$mqtl2)
#df$mqtl2 <- factor(df$mqtl, levels = c("Extreme Height SNPs (n=60)","Extreme Height mQTL (n=48)","Extreme Height mQTL SDS (n=4)","Height SNPs (n=3290)","Height mQTL(n=3948)","Height mQTL SDS (n=44)","Schizophrenia SNPs (n=159)","Schizophrenia mQTL (n=133)","Schizophrenia mQTL SDS (n=9)","CVD SNPs (n=202)","CVD mQTL (n=15)","CVD mQTL SDS (n=2)"))

p1<-ggplot(df, aes(y=es, x=mqtl2)) +
geom_jitter(shape=16, position=position_jitter(0.2)) +
facet_wrap(~ trait,scales="free",nrow=1) +
labs(y="Genetic Variance",x="mQTL category") +

stat_summary(
  aes(y=es, x=mqtl2,color="red"),
  fun.data="mean_sdl",  fun.args = list(mult=1), 
  geom = "pointrange",  size = 0.4)+
theme(legend.position = 'none') +
theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(p1,file="FigSX_GeneticVariance_boxplot.pdf",height=6,width=10)

load("/newshared/godmc/database_files/snps.rdata")
df2<-df[which(df$SNP%in%out_df2$name),]

nrow(df)
#[1] 9001
nrow(df2)
#[1] 8943

table(df2$mqtl2)
df2$mqtl2<-gsub("All SNPs (n=197)","All SNPs (n=185)",df2$mqtl2,fixed=T)
df2$mqtl2<-gsub("All SNPs (n=3229)","All SNPs (n=3210)",df2$mqtl2,fixed=T)
df2$mqtl2<-gsub("All SNPs (n=80)","All SNPs (n=79)",df2$mqtl2,fixed=T)
df2$mqtl2<-gsub("mQTL (n=4910)","mQTL (n=4888)",df2$mqtl2,fixed=T)
df2$mqtl2<-gsub("mQTL (n=55)","mQTL (n=54)",df2$mqtl2,fixed=T)
df2$mqtl2<-gsub("mQTL (n=75)","mQTL (n=73)",df2$mqtl2,fixed=T)
df2$mqtl2<-gsub("mQTL (n=152)","mQTL (n=151)",df2$mqtl2,fixed=T)


p1<-ggplot(df2, aes(y=es, x=mqtl2)) +
geom_jitter(shape=16, position=position_jitter(0.2)) +
facet_wrap(~ trait,scales="free",nrow=1) +
labs(y="Genetic Variance",x="mQTL category") +

stat_summary(
  aes(y=es, x=mqtl2,color="red"),
  fun.data="mean_sdl",  fun.args = list(mult=1), 
  geom = "pointrange",  size = 0.4)+
theme(legend.position = 'none') +
theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(p1,file="FigSX_GeneticVariance_boxplot_filtered.pdf",height=6,width=10)




#p1<-ggplot(df, aes(y=es, x=mqtl2)) +
#  geom_dotplot(binaxis = "y", stackdir = "center",fill = "white") +
#  facet_wrap(~ trait,scales="free",nrow=1) +
#  labs(y="Genetic Variance",x="mQTL category") +
#  theme(axis.text.x=element_text(angle=90,hjust=1))
#ggsave(p1,file="FigSX_GeneticVariance_boxplot.pdf",height=6,width=10)

 