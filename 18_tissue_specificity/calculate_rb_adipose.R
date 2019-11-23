load("tissue.RData")
load("~/repo/godmc_phase2_analysis/07_enrichments/mean_allcpgs.Robj")
df.all$meancat<-NA
hypo<-which(df.all$meancpg<0.2)
df.all$meancat[hypo]<-"hypo"
intermediate<-which(df.all$meancpg>0.2&df.all$meancpg<0.8)
df.all$meancat[intermediate]<-"intermediate"
hyper<-which(df.all$meancpg>0.8)
df.all$meancat[hyper]<-"hyper"

clumped<-merge(clumped,df.all[,c("cpg","meancpg","meancat")],by.x="cpg",by.y="cpg",all.x=T)

 #standarized beta and SE in SD unit
# z represent z statistics, p presents allele frequency, n represents sample size.\

#
#            r        se_r
#[1,] 0.880942 0.008752756


#Qi et al:The mean brain–blood r b estimate from two samples was 0.78 (s.e. = 0.006) 

calcu_std_b_se<-function(z,p,n){
    std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
    std_se=1/sqrt(2*p*(1-p)*(n+z^2))
    res<-data.frame(std_b_hat,std_se);
    return(res)
}

# calcualte rb
# b1 and se1 represent the estimate and SE for top eQTLs across probes in one tissue, b2 and se2 represent the estimate and SE for top eQTLs across probes in the other tissue
# theta = sample overlap * phnotypic correlation; theta could also be estimated from null SNPs; if two samples were independent, theta = 0.
# Please note that the effect allele of one SNP between two tissues should be the same.
calcu_cor_true<-function(b1,se1,b2,se2,theta){
    idx=which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
    if(length(idx)>0){
        b1=b1[-idx];se1=se1[-idx]
        b2=b2[-idx];se2=se2[-idx]
        theta=theta[-idx]
    }

    var_b1=var(b1,na.rm=T)-mean(se1^2,na.rm=T)
    var_b2=var(b2,na.rm=T)-mean(se2^2,na.rm=T)
    cov_b1_b2=cov(b1,b2,use="complete.obs")-mean(theta,na.rm=T)*sqrt(mean(se1^2,na.rm=T)*mean(se2^2,na.rm=T))
    r=cov_b1_b2/sqrt(var_b1*var_b2)

    r_jack=c()
    n=length(b1)
    for(k in 1:n) {
       b1_jack=b1[-k];se1_jack=se1[-k];var_b1_jack=var(b1_jack,na.rm=T)-mean(se1_jack^2,na.rm=T)
       b2_jack=b2[-k];se2_jack=se2[-k];var_b2_jack=var(b2_jack,na.rm=T)-mean(se2_jack^2,na.rm=T)
       theta_jack=theta[-k];
       cov_e1_jack_e2_jack=mean(theta_jack,na.rm=T)*sqrt(mean(se1_jack^2,na.rm=T)*mean(se2_jack^2,na.rm=T))
       cov_b1_b2_jack=cov(b1_jack,b2_jack,use="complete.obs")-cov_e1_jack_e2_jack
       r_tmp=cov_b1_b2_jack/sqrt(var_b1_jack*var_b2_jack)
       r_jack=c(r_jack,r_tmp)
    }
    r_mean=mean(r_jack,na.rm=T)
    idx=which(is.na(r_jack))
    if(length(idx)>0){
        se_r=sqrt((n-1)/n*sum((r_jack[-idx]-r_mean)^2))
    }else{
      se_r=sqrt((n-1)/n*sum((r_jack-r_mean)^2))
    }

    res<-cbind(r,se_r)
    return(res)
}

m<-match(a2$MARKERNAME,clumped$id,)
cl2<-clumped[m,]

cat<-unique(cl2$cpg_cis2)

cor.out<-data.frame()
for (i in 1:length(cat)){
cat(cat[i],"\n")
w<-which(cl2$cpg_cis2==cat[i])
out<-calcu_cor_true(b1=cl2$Effect[w],se1=cl2$StdErr[w],b2=a2$BETA2[w],se2=a2$SE[w],theta=0)
cor.out<-rbind(cor.out,unlist(out))
}
cor.out<-data.frame(cor.out,category=cat)
save(cor.out,file="rb_adipose.Robj")

###
w<-which(cl2$pval<1e-14)
cl2<-cl2[w,]
a2<-a2[w,]

cat<-unique(cl2$cpg_cis2)
cor.out<-data.frame()
for (i in 1:length(cat)){
cat(cat[i],"\n")
w<-which(cl2$cpg_cis2==cat[i])
out<-calcu_cor_true(b1=cl2$Effect[w],se1=cl2$StdErr[w],b2=a2$BETA2[w],se2=a2$SE[w],theta=0)
cor.out<-rbind(cor.out,unlist(out))
}
cor.out<-data.frame(cor.out,category=cat)
save(cor.out,file="rb.adipose.filtered.Robj")
###
hla<-which(cl2$snpchr=="chr6"&cl2$snppos>29570005&cl2$snppos<33377657)
cl2<-cl2[-hla,]
a2<-a2[-hla,]

cat<-unique(cl2$cpg_cis2)
cor.out<-data.frame()
for (i in 1:length(cat)){
cat(cat[i],"\n")
w<-which(cl2$cpg_cis2==cat[i])
cl3<-cl2[w,]

####if(cat[i]=="cisonly"|cat[i]=="cis+trans_cis"){
o<-order(cl3$pval)
cl3<-cl3[o,]
m<-match(unique(cl3$cpg),cl3$cpg)
cl3<-cl3[m,]
m<-match(cl3$id,a2$MARKERNAME)
a3<-a2[m,]

####}

out<-calcu_cor_true(b1=cl3$Effect,se1=cl3$StdErr,b2=a3$BETA2,se2=a3$SE,theta=0)
cor.out<-rbind(cor.out,unlist(out))
}
cor.out<-data.frame(cor.out,category=cat)
save(cor.out,file="rb.adipose.filtered_primary.Robj")
###
cat2<-unique(cl2$cis)
cor.out<-data.frame()
for (i in 1:length(cat2)){
cat(cat2[i],"\n")
w<-which(cl2$cis==cat2[i])
cl3<-cl2[w,]

o<-order(cl3$pval)
cl3<-cl3[o,]
m<-match(unique(cl3$cpg),cl3$cpg)
cl3<-cl3[m,]
m<-match(cl3$id,a2$MARKERNAME)
a3<-a2[m,]

out<-calcu_cor_true(b1=cl3$Effect,se1=cl3$StdErr,b2=a3$BETA2,se2=a3$SE,theta=0)
cor.out<-rbind(cor.out,unlist(out))
}
cor.out<-data.frame(cor.out,category=cat2)
save(cor.out,file="rb.adipose.filtered_primary_cistrans.Robj")
###

cl2$cpg_cis3<-paste0(cl2$cpg_cis2,"_",cl2$meancat)

cat<-unique(cl2$cpg_cis3)

cor.out<-data.frame()
for (i in 1:length(cat)){
cat(cat[i],"\n")
w<-which(cl2$cpg_cis3==cat[i])
cl3<-cl2[w,]

####if(cat[i]=="cisonly"|cat[i]=="cis+trans_cis"){
o<-order(cl3$pval)
cl3<-cl3[o,]
m<-match(unique(cl3$cpg),cl3$cpg)
cl3<-cl3[m,]
m<-match(cl3$id,a2$MARKERNAME)
a3<-a2[m,]

####}

out<-calcu_cor_true(b1=cl3$Effect,se1=cl3$StdErr,b2=a3$BETA2,se2=a3$SE,theta=0)
cor.out<-rbind(cor.out,unlist(out))
}
cor.out<-data.frame(cor.out,category=cat)
save(cor.out,file="rb.adipose.filtered_primary_hypo.Robj")

gc()



#> out
#             r       se_r
#[1,] 0.6225506 0.01612748
#> cat[i]
#[1] "cis+trans_trans"


