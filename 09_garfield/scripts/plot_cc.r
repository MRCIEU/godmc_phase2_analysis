library(ggplot2)
library(gridExtra)
library(grid)

path="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/output"
l<-list.files(path,pattern="27863252")
meta<-read.table("/projects/MRC-IEU/research/projects/ieu2/p5/021/working/data/Astle2016/meta_data_forMRBase_Astle2016.txt",he=T,sep="\t")
spl<-strsplit(as.character(l),split="-")
spl<-do.call("rbind",spl)
m<-match(spl[,2],meta$Study.accession)
meta<-meta[m,]

cc<-data.frame()
for (i in 1:length(l)){

r<-read.table(paste0(path,"/",l[i],"/garfield.test.",l[i],".out"),he=T)
cat(nrow(r),"\n")
r$Celltype<-rep(meta$Reported.trait[i],4)
r$Type<-rep(meta$Study.accession[i],4)
cc<-rbind(cc,r)
}

padj<-p.adjust(cc$Pvalue,method="BH")
cc$padj<-padj

cc$Annotation<-gsub("cis","cis only",cc$Annotation)
cc$Annotation<-gsub("trans","trans only",cc$Annotation)
cc$Annotation<-gsub("cis onlytrans only","cis+trans",cc$Annotation)

p1<-ggplot(cc, aes(fill=Annotation, y=-log10(cc$Pvalue), x=Celltype)) + 
    geom_bar(position="dodge", stat="identity") +
    labs(x="",y="-log10 (Pvalue)") +
    theme(axis.text.x = element_blank())
    scale_fill_brewer(type="qual")

m<-max(cc$Beta)+1

p2<-ggplot(cc, aes(color=Annotation, y=Beta, x=Celltype)) +
 scale_fill_brewer(type="qual") +
    #geom_bar(position="dodge", stat="identity") +
    geom_point(position=position_dodge(.9)) +
    geom_errorbar(aes(x=Celltype,ymin=CI95_lower, ymax=CI95_upper), width=.2,position=position_dodge(.9)) +
    #scale_fill_brewer(type="qual") +
    labs(x="Phenotype",y="log OR") +
    #scale_y_continuous(trans = 'log10',breaks=c(.2,.3,.4,.5,1:m),limits=c(0.2,m)) +
    scale_y_continuous(limits=c(-2,6)) +
    geom_hline(aes(yintercept = 0)) +
    #scale_fill_brewer(type="qual") +
    theme(axis.text.x=element_text(angle = 90, hjust = 10))

#ggsave(p1,file="cellcounts_plot.pdf",height=6,width=10)

pdf("GWA_enrichment_cellcounts.pdf",height=10,width=18)
#create layout, assign it to viewport, push viewport to plotting device
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2:3, 1))

dev.off()
