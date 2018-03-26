library(GenomicRanges)
library(data.table)
library(tidyr)
arguments<-commandArgs(T)
i<-as.numeric(arguments[1])


path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/encode_tfbs"
l<-list.files(path)
chr<-paste0("chr",i)

r<-read.table(paste0("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/maftssd/",chr))
head(r)
r<-data.frame(snpchr=chr,min=as.numeric(r$V1),max=as.numeric(r$V1),r)
#       V1           V2      V3 V4        V5          V6
#1 9411245 0.0017005394 -271945  0 0.2287712 0.003623188
#2 9411318 0.1636909693 -271872  0 0.2287712 0.003623188
#3 9411354 0.0017026922 -271836  0 0.2287712 0.003623188
#4 9411370 0.0001440761 -271820  0 0.2287712 0.003623188
#5 9411384 0.0013298947 -271806  0 0.2287712 0.003623188
#6 9411540 0.0008329048 -271650  0 0.2287712 0.003623188

#collapse overlaps
gr_snp = with(r,GRanges(seqnames=snpchr,ranges=IRanges(min,max),strand=Rle("+")))

df<-data.frame(r$min)

ov<-r$min
for (i in 1:length(l)){
#for (i in 1:10){
cat(i,"\n")
bed<-read.table(paste0(path,"/",l[i]))
bed$V1<-gsub("chrX","chr23",bed$V1)
bed<-bed[which(bed$V1==chr),]
#     V1    V2    V3   V4    V5 V6
#1 chr21 11529 11731 CTCF Dnd41  *
#2 chr21 11619 11929 CTCF Dnd41  *
#3 chr21 11668 11978 CTCF Dnd41  *
#4 chr21 13183 13396 CTCF Dnd41  *
#5 chr21 14423 14580 CTCF Dnd41  *
#6 chr21 18733 18929 CTCF Dnd41  *

if(nrow(bed)>0){
bed$V2<-as.numeric(bed$V2)
bed$V3<-as.numeric(bed$V3)
bed<-unique(data.frame(chr=bed$V1,start=bed$V2,end=bed$V3,strand="*",antibody=bed$V4,celltype=bed$V5))
gr_range<-makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE,starts.in.df.are.0based=TRUE) 
overlap<-countOverlaps(gr_snp,gr_range)
print(length(which(overlap==1)))
w<-which(overlap>0)
if(length(w)>0){overlap[w]<-1}
}

if(nrow(bed)==0){overlap<-rep(0,nrow(r))}

ov<-data.frame(ov,overlap)
}
cols <- names(ov)[-1]
#ov$ann <- apply( ov[ , cols ] , 1 , paste , collapse = "" )
#ov <- ov[ , !(names(ov)%in%cols) ]

#ov$ann <- do.call(paste, c(ov[cols], sep=""))

ov2<-unite(ov, ann, -1,sep="") 

a<-apply(ov2,1,function(x) nchar(x[2]))
print(length(table(a)))



write.table(ov2,paste0("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/annotation-tfbs/",chr),sep=" ",quote=FALSE,col.names=F,row.names=F)

#0 or 1-based?
#test<-which(r$min%in%bed$V2)
#test<-r[test[1],"min"]
#w1<-which(bed$V2<=test)
#w2<-which(bed$V3>test)
#w<-which(w1%in%w2)

#which(df$r.min%in%test)
#df[8187,]
#        r.min overlap
#8187 14926178       2

#16051237 00000
#16051249 00000
#16051477 00000

#Index Annotation Celltype Tissue Type Category
#0 ihs blood blood mqtls mqtl
#1 fst blood blood mqtls mqtl
#2 xpehhchb blood blood mqtls mqtl
#3 xpehhyri blood blood mqtls mqtl
#4 sds blood blood mqtls mqtl


#R CMD BATCH --no-save --no-restore '--args '$i'' preparegarfield_selection_sds_trans.R preparegarfield_selection${i}_sds_trans.Rout
