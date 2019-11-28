library(dplyr)
df<-read.table("adipose_rb.txt",sep="\t",he=T)
r<-read.table("brain_rb.txt",sep="\t",he=T)
df<-merge(df[,c("category","n","Rb")],r[,c("category","n","Rb")],by.x="category",by.y="category")
names(df)<-c("Category","N assoc Adipose","Rb Adipose","N assoc brain","Rb Brain")

load("tissue.RData")

gp1<-clumped %>% group_by(cis) %>% dplyr::summarize(count = n_distinct(cpg))
names(gp1)[1]<-"Category"
gp1$Category<-gsub("FALSE","trans any",gp1$Category)
gp1$Category<-gsub("TRUE","cis any",gp1$Category)

gp2<-clumped %>% group_by(cpg_cis2) %>% dplyr::summarize(count = n_distinct(cpg))
names(gp2)[1]<-"Category"

gp3<-clumped %>% group_by(cpg_cis3) %>% dplyr::summarize(count = n_distinct(cpg))
names(gp3)[1]<-"Category"

blood<-rbind(gp1,gp2,gp3)

names(blood)[2]<-"N assoc blood"

df<-merge(blood,df,by.x="Category",by.y="Category")

df$adipose_blood_ratio<-round(df[,3]/df[,2],2)
df$brain_blood_ratio<-round(df[,5]/df[,2],2)

write.table(df,"adipose_brain_rb.txt",sep="\t",col.names=T,row.names=F,quote=F)

gp3<-clumped %>% group_by(cpg_cis3) %>% dplyr::summarize(frq = mean(Freq1), I2=mean(HetISq),Effect=mean(abs(Effect)))

cl = group_by(clumped,cpg_cis3) %>% dplyr::summarize(count = n_distinct(cpg))

cat<-unique(clumped$cpg_cis3)

cl.new<-data.frame()
for (i in 1:length(cat)){
cat(cat[i],"\n")
w<-which(clumped$cpg_cis3==cat[i])
cl2<-clumped[w,]
o<-order(cl2$pval)
cl2<-cl2[o,]
m<-match(unique(cl2$cpg),cl2$cpg)
cl2<-cl2[m,]
cl.new<-rbind(cl.new,cl2)}

gp3.new<-cl.new %>% group_by(cpg_cis3) %>% dplyr::summarize(frq = mean(Freq1), I2=mean(HetISq),Effect=mean(abs(Effect)))
very low = promoter, low-intermediate = enhancer, high = quiescent