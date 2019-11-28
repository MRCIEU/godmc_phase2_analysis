#module add languages/R-3.5.1-ATLAS-gcc-6.1
#mkdir /panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation
#mkdir /panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation/7regions
#mkdir /panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation/7regions/regions
path="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation/7regions"
#tss<-read.table("/newshared/godmc/gencode/gencode.v27lift37.annotation.gtf.gz",skip=5,sep="\t")
#names(tss)[1:5]<-c("chromosome","source","type","start","end")

library(GenomicRanges)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(annotatr)
library(data.table)


Illumina450=IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
Illumina450 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Illumina450_dt=as.data.table(Illumina450[,1:3])
Illumina450_dt[,cpgID:=row.names(Illumina450),]
Illumina450_dt[,cpgstart_pre:=ifelse(strand=="-",pos-500,pos-499),]
Illumina450_dt[,cpgend_pre:=ifelse(strand=="-",pos+500,pos+501),]

b<-builtin_annotations()

annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
    'hg19_genes_intronexonboundaries')
annotations = build_annotations(genome = 'hg19', annotations = annots)

table(annotations$type)

#                 hg19_cpg_inter                hg19_cpg_islands 
#                          20613                           28691 
#               hg19_cpg_shelves                 hg19_cpg_shores 
#                          43747                           51914 
#              hg19_genes_1to5kb                hg19_genes_3UTRs 
#                          82960                           64362 
#               hg19_genes_5UTRs                hg19_genes_exons 
#                         114305                          742493 
#          hg19_genes_intergenic hg19_genes_intronexonboundaries 
#                          17027                          690906 
#             hg19_genes_introns            hg19_genes_promoters 
#                         659327                           82960 


#3'UTR, 5'UTR, Downstream, Exon, Intergenic, Intron, Upstream

df<-data.frame(annotations)
#w<-which(df$strand=="-")
w<-which(df$end<df$start)
start<-df$start
start[w]<-df$end[w]
stop<-df$end
stop[w]<-df$start[w]

res<-unique(data.frame(chr=df$seqnames,start,stop,type=df$type))
type<-unique(df$type)

meta.out<-data.frame()
for (i in 1:length(type)){
res2<-res[which(res$type==type[i]),]
file.name<-paste0(type[i],".bed")
meta<-data.frame(filename=file.name,description=type[i],tissue=NA,antibody=NA,treatment=NA,dataSource="hg19_genes")
write.table(res2,paste(path,"regions",file.name,sep="/"),sep="\t",quote=F,row.names=F,col.names=F)
meta.out<-rbind(meta.out,meta)
}

write.table(meta.out,paste0(path,"/index.txt"),sep="\t",col.names=T,row.names=F,quote=F)

coll<-data.frame(collector="Josine Min",date="2019-09-20",source="AnnotationHub_2.14.5",description="hg19_genes")
write.table(coll,paste0(path,"/collection.txt"),sep="\t",col.names=T,row.names=F,quote=F)


db=loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation/")

