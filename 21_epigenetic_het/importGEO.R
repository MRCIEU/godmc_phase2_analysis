source("~/LUMCsamples/geo.r")
library(data.table)
library(meffil)
library(lumi)
library(gplots)
library(gridGraphics)
library(grid)
library(gridExtra)
library(ggplot2)

#original gplots package is very slow due to a bug
#remove.packages('gplots')
#library('devtools')
#install_github("ChristophH/gplots")

## dataset monocytes
samples <- geo.samples("GSE56046")
samples <- with(samples, {
    data.frame(gsm=geo_accession,
               gse="GSE56046",
               age=get.characteristic(characteristics_ch1,"age"),
               racegendersite=get.characteristic(characteristics_ch1,"racegendersite"),
               ChIP=get.characteristic(characteristics_ch1,"ChIP"),
               well=get.characteristic(characteristics_ch1,"well"),
               bcell=get.characteristic(characteristics_ch1,"bcell"),
               tcell=get.characteristic(characteristics_ch1,"tcell"),
               nkcell=get.characteristic(characteristics_ch1,"nkcell"),
               neutro=get.characteristic(characteristics_ch1,"neutro"),
               plaque=get.characteristic(characteristics_ch1,"plaque"),
               cac=get.characteristic(characteristics_ch1,"plaque"),
               stringsAsFactors=F)
})
#d<-fread("zcat ~/MESA/GSE56046.tmp.gz") #485578 #doesn't work properly



d<-fread("~/MESA/GSE56046.tmp")
d[1:5,1:5]
#       ID_REF 100001.Mvalue 100001.detectionPval 100002.Mvalue
#1: cg00000029    -0.2482799                    0    -0.0203916
#2: cg00000108     3.5163442                    0     3.8562861
#3: cg00000109     2.7526248                    0     3.2674936
#4: cg00000165    -2.0308491                    0    -1.4621109
#5: cg00000236     1.0695492                    0     1.1464122
#   100002.detectionPval
#1:                    0
#2:                    0
#3:                    0
#4:                    0
#5:                    0

g<-grep("Mvalue",colnames(d))
d<-data.frame(d)
dm<-d[,g]
db<-m2beta(dm)
colnames(db)<-gsub("Mvalue","mono",colnames(db))
row.names(db)<-d[,1]

message("Identifying methylation outliers")

	niter <- 3
	outlier_threshold <- 10
	db.copy <- db
	for(i in 1:niter)
	{
		sds <- rowSds(db.copy, na.rm=T)
		means <- rowMeans(db.copy, na.rm=T)
		db.copy[db.copy > means + sds*outlier_threshold | db.copy < means - sds*outlier_threshold] <- NA
	}

	db[db.copy] <- NA
	outlier_count <- apply(db.copy, 1, function(x) sum(is.na(x)))
	table(outlier_count)
	db.copy <- is.na(db.copy)
	db[db.copy] <- NA
	outlier_count <- apply(db, 1, function(x) sum(is.na(x)))
	table(outlier_count)

rowm<-data.frame(cpg=d[,1],monocytes=rowMeans(db))

#tcell<-fread("zcat ~/MESA/GSE56581.tmp.gz") #doesn't work properly
tcell<-fread("~/MESA/GSE56581.tmp")

g<-grep("Mvalue",colnames(tcell))
tcell<-data.frame(tcell)
dt<-tcell[,g]
tb<-m2beta(dt)
colnames(tb)<-gsub("Mvalue","tcell",colnames(tb))
row.names(tb)<-tcell[,1]

message("Identifying methylation outliers")

	niter <- 3
	outlier_threshold <- 10
	tb.copy <- tb
	for(i in 1:niter)
	{
		sds <- rowSds(tb.copy, na.rm=T)
		means <- rowMeans(tb.copy, na.rm=T)
		tb.copy[tb.copy > means + sds*outlier_threshold | tb.copy < means - sds*outlier_threshold] <- NA
	}
	
	tb[tb.copy] <- NA
	outlier_count <- apply(tb.copy, 1, function(x) sum(is.na(x)))
	table(outlier_count)
	tb.copy <- is.na(tb.copy)
	db[db.copy] <- NA
	outlier_count <- apply(db, 1, function(x) sum(is.na(x)))
	table(outlier_count)

	#index <- which(is.na(tb), arr.ind = TRUE)
	#if (length(index)>0){
    #message("Replace ",length(index)," missing values with rowmeans")
    #tb[index] <- rowMeans(tb, na.rm = TRUE)[index[, "row"]] }



rowm2<-data.frame(cpg=tcell[,1],Tcells=rowMeans(tb))



##
load("../07_enrichments/mean_allcpgs.Robj")

ind_m<-list()
ind_t<-list()
ct_list<-list()

cat<-unique(df.all$cpg_cis)
for (i in 1:length(cat)){
ct<-df.all[which(df.all$cpg_cis==cat[i]),c("cpg","meancpg")]
m<-match(ct$cpg,rowm$cpg)
m2<-match(ct$cpg,rowm2$cpg)
ct2<-data.frame(blood=ct[,2],monocytes=rowm[m,c("monocytes")],Tcells=rowm2[m2,c("Tcells")])
row.names(ct2)<-ct$cpg
ct_ind_m<-data.frame(db[m,])
ct_ind_t<-data.frame(tb[m2,])
row.names(ct_ind_m)<-ct$cpg
row.names(ct_ind_t)<-ct$cpg

o<-order(ct2$blood,decreasing=TRUE)
ct2<-ct2[o,]
ord<-row.names(ct2)[o]
ct2<-as.matrix(ct2)

ct_ind_m<-ct_ind_m[o,]
ct_ind_m<-as.matrix(ct_ind_m)
ct_ind_t<-ct_ind_t[o,]
ct_ind_t<-as.matrix(ct_ind_t)

ind_t[[i]]<-ct_ind_t
ind_m[[i]]<-ct_ind_m
ct_list[[i]]<-ct2
}

#save(ct_ind_m,ct_ind_t,ct2,cat,file="heatmap_monocytes_tcells.Robj")
save(ind_m,ind_t,ct_list,cat,file="heatmap_monocytes_tcells.Robj")
###

col = c("red","orange","yellow")
breaks <- c(0, 0.20, 0.80, 1)

for (i in 1:length(cat)){
ct2<-as.matrix(ct_list[[i]])
pdf(paste0("heatmap.",cat[i],".pdf"),height=6,width=6)
heatmap.2(ct2,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=FALSE,Rowv=FALSE,labRow=FALSE,labCol=colnames(ct2),breaks=breaks,col=col)
dev.off()

ct2<-data.frame(ct2)
length(which(ct2$blood>0.2&ct2$blood<0.8))
length(which(ct2$blood>0.2&ct2$blood<0.8))/nrow(ct2) #66% 7847
length(which(ct2$monocytes>0.2&ct2$monocytes<0.8))/nrow(ct2) #62%
length(which(ct2$Tcells>0.2&ct2$Tcells<0.8)) #7686 65%

length(which(ct2$blood>0.2&ct2$blood<0.8&ct2$monocytes>0.2&ct2$monocytes<0.8& ct2$Tcells>0.2&ct2$Tcells<0.8))/nrow(ct2)
length(which(ct2$blood<0.2))/nrow(ct2) #2599 22%
length(which(ct2$blood>0.8))/nrow(ct2) #1456 12%

pdf(paste0("heatmap_ind",cat[i],".pdf"),height=6,width=8)
heatmap.2(ind_t[[i]],keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="column",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=T)
heatmap.2(ind_m[[i]],keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="column",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=T)
dev.off()

}

#p1<-heatmap.2(ct2,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=FALSE,Rowv=FALSE,labRow=F,labCol=F)
#p2<-heatmap.2(ct_ind_t,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=F)
#p3<-heatmap.2(ct_ind_m,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=F)
#l<-list(p1,p2,p3)

#ggMultiPdf <- function(plot_list, filename){
 
#plots_gr <- marrangeGrob(plot_list, nrow = 1, ncol = 3)
 
#ggsave(filename = filename,
#       plot = plots_gr,
#       dpi = 60,
#       device = cairo_pdf,
#       onefile = TRUE)
 
#}
 
#ggMultiPdf(l,"heatmap_ind.pdf")

#pdf(file="heatmap_ind.pdf",height=6,width=8)


#library(gridGraphics)
#grab_grob <- function(){
#  grid.echo()
#  grid.grab()
#}

#g <- grab_grob()
#grid.newpage()

#lay <- grid.layout(nrow = 1, ncol=3)
#pushViewport(viewport(layout = lay))

#pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
#heatmap.2(ct2,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=FALSE,Rowv=FALSE,labRow=F,labCol=F)
#upViewport()

#pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
#heatmap.2(ct_ind_t,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=F)
#upViewport()

#pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
#heatmap.2(ct_ind_m,keysize=1.5,key.title=NA,key.xlab="mean DNAm",dendrogram="none",trace="none",cexCol=1,Colv=TRUE,Rowv=FALSE,labRow=F,labCol=F)
#upViewport()

#upViewport()
#dev.off()

