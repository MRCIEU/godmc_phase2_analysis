library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stats)
library(missMethyl)
library(KEGGREST)

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
allcpgs <- rownames(ann)

cpgs<-read.table("../results/intracranial_volume.txt",sep="\t",he=F,colClasses="character")
go.res <- gometh(sig.cpg=cpgs[,1], all.cpg=allcpgs,array.type=c("450K"))
go.res[go.res$FDR<0.05,]

go.res2 <- gometh(sig.cpg=cpgs[,1], all.cpg=allcpgs,array.type=c("450K"),collection="KEGG")
go.res2[which(go.res2$FDR<0.05),]

save(go.res,go.res2,file="../results/goenrichments_intracranialvolume.rdata")


cpgs<-read.table("../results/intracranial_volume2.txt",sep="\t",he=F,colClasses="character")
go.res <- gometh(sig.cpg=cpgs[,1], all.cpg=allcpgs,array.type=c("450K"))
go.res[go.res$FDR<0.05,]

go.res2 <- gometh(sig.cpg=cpgs[,1], all.cpg=allcpgs,array.type=c("450K"),collection="KEGG")
go.res2[which(go.res2$FDR<0.05),]

save(go.res,go.res2,file="../results/goenrichments_intracranialvolume2.rdata")
