library(genefilter)
library(gplots)

load("norm.beta.pc15.Robj")

message("Identifying methylation outliers")

	niter <- 3
	outlier_threshold <- 10
	norm.beta.copy <- norm.beta
	for(i in 1:niter)
	{
		sds <- rowSds(norm.beta.copy, na.rm=T)
		means <- rowMeans(norm.beta.copy, na.rm=T)
		norm.beta.copy[norm.beta.copy > means + sds*outlier_threshold | norm.beta.copy < means - sds*outlier_threshold] <- NA
	}

	norm.beta[norm.beta.copy] <- NA
	outlier_count <- apply(norm.beta.copy, 1, function(x) sum(is.na(x)))
	table(outlier_count)
	norm.beta.copy <- is.na(norm.beta.copy)
	norm.beta[norm.beta.copy] <- NA
	outlier_count <- apply(norm.beta, 1, function(x) sum(is.na(x)))
	table(outlier_count)

samples<-read.table("~/LUMCsamples/samplesheet.GEO.txt",sep="\t",he=T)
names(samples)<-gsub("sex","Sex",names(samples))
names(samples)<-gsub("sentrix_id","Sample_Name",names(samples))
slide<-do.call("rbind",strsplit(as.character(samples$Sample_Name),split="_"))
samples<-data.frame(samples,slide=slide[,1])
samples$tissue<-gsub(".Buffy","",samples$tissue)
samples$tissue<-gsub(".Skel","",samples$tissue)
samples$tissue<-gsub("Left.","",samples$tissue)
samples$tissue<-gsub(".Sub","",samples$tissue)

m<-match(colnames(norm.beta),samples$Sample_Name,)
samples<-samples[m,]

rowt<-data.frame(cpg=row.names(norm.beta))
cat<-unique(samples$tissue)
for (i in 1:length(cat)){
cat(i,"\n")
w<-which(samples$tissue%in%cat[i])
tiss<-norm.beta[,w]
tiss<-data.frame(tiss=rowMeans(tiss))
names(tiss)[1]<-cat[i]
rowt<-data.frame(rowt,tiss)
}

###


load("../07_enrichments/mean_allcpgs.Robj")

ind_m<-list()
ind_t<-list()
ct_list<-list()

cat<-unique(df.all$cpg_cis)
for (i in 1:length(cat)){
ct<-df.all[which(df.all$cpg_cis==cat[i]),c("cpg","meancpg")]
m<-match(ct$cpg,rowt$cpg)
ct2<-data.frame(blood=ct[,2],rowt[m,-1])
row.names(ct2)<-ct$cpg
ct_ind_t<-data.frame(norm.beta[m,])
row.names(ct_ind_t)<-ct$cpg

o<-order(ct2$blood,decreasing=TRUE)
ct2<-ct2[o,]
ord<-row.names(ct2)[o]
ct2<-as.matrix(ct2)

ct_ind_t<-ct_ind_t[o,]
ct_ind_t<-as.matrix(ct_ind_t)

ind_t[[i]]<-ct_ind_t
ct_list[[i]]<-ct2
}

#save(ct_ind_m,ct_ind_t,ct2,cat,file="heatmap_monocytes_tcells.Robj")
save(ind_t,ct_list,cat,file="tissue.Robj")
###

col = c("red","orange","yellow")
breaks <- c(0, 0.20, 0.80, 1)

library(meffil)
y<-meffil.get.features("450k")
y<-data.frame(cpg=y$name,cpgchr=y$chromosome)

w<-which(cat%in%"cis only")
data<-data.frame(ct_list[[w]])
m<-match(row.names(data),y$cpg)
y2<-y[m,]
hyper<-which(data$blood>0.8)
table(y2[hyper,2])

hyper<-matrix(nr=4,nc=4)
for (i in 1:length(cat)){
ct2<-as.matrix(ct_list[[i]])
colnames(ct2)[1]<-"Blood"
colnames(ct2)[2]<-"PBMC"


png(paste0("heatmap.tissue.",cat[i],".png"),height=600,width=600)
heatmap.2(ct2,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=FALSE,Rowv=FALSE,labRow=FALSE,labCol=colnames(ct2),breaks=breaks,col=col)
dev.off()

ct2<-data.frame(ct2)
hyper[i,1]<-nrow(ct2)

hyper[i,3]<-length(which(ct2$Blood>0.2&ct2$Blood<0.8))
#length(which(ct2$Blood>0.2&ct2$Blood<0.8))/nrow(ct2))#66% 7847
hyper[i,2]<-length(which(ct2$Blood<0.2)) #2599 22%
hyper[i,4]<-length(which(ct2$Blood>0.8)) #1456 12%

#pdf(paste0("heatmap.tissue_ind.",cat[i],".pdf"),height=6,width=8)
#heatmap.2(ind_t[[i]],keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="column",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=colnames(ct2))
#dev.off()

}

colnames(hyper)<-c("N","hypo","intermediate","hyper")
hyper<-data.frame(category=cat,hyper)

for (i in 1:length(cat)){
ct2<-as.matrix(ct_list[[i]])
colnames(ct2)[1]<-"Blood"
colnames(ct2)[2]<-"PBMC"

load("~/repo/godmc_phase2_analysis/07_enrichments/mean_allcpgs.Robj")
df.all$meancat<-NA
hypo<-which(df.all$meancpg<0.2)
df.all$meancat[hypo]<-"hypo"
intermediate<-which(df.all$meancpg>0.2&df.all$meancpg<0.8)
df.all$meancat[intermediate]<-"intermediate"
hyper<-which(df.all$meancpg>0.8)
df.all$meancat[hyper]<-"hyper"

m<-match(row.names(ct2),df.all$cpg)
ct2<-data.frame(ct2,meancat=NA)
ct2$meancat<-df.all$meancat[m]
ct3<-ct2[,1:(ncol(ct2)-1)]
pdf(paste0("meanDNAm.tissue.",cat[i],".pdf"),height=6,width=6)

col = c("red","orange","yellow")
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, cex=0.5,col = col[ct2$meancat])
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
}
# Create the plots
pairs(ct3, lower.panel = NULL, pch = 19, cex=0.5,col = col[ct2$meancat])

dev.off()

