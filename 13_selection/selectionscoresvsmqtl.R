load("../results/enrichments/snpcontrolsets_selection.rdata")
f.all1<-f.all
load("../results/enrichments/snpcontrolsets.rdata")
df<-data.frame(f.all,f.all1[,c(29:43)])


library(ggplot2)

p1<-ggplot(df, aes(sds_score, max_abs_Effect)) +
geom_point() +
theme_minimal()
ggsave(p1, file="sds.pdf", width=10, height=10)

df$min_pval2<-df$min_pval
w<-which(df$min_pval==0)
df$min_pval2[w]<-min(df$min_pval[-w],na.rm=T)

p2<-ggplot(df, aes(sds_score, -log10(df$min_pval2))) +
geom_point() +
theme_minimal()
ggsave(p2, file="sdspval.pdf", width=10, height=10)