setwd("/data/home/eptgr/GoDMC/")
library("cmapR")
library("data.table")
library("dplyr")
library("biomaRt")

#source("convertIDs.R")

## get indices of relevant perturbations
ds <- parse.gctx("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", rid=1)
cids <- as.data.frame(ds@cid)
colnames(cids) <- c('sig_name')
cids$index <- row.names(cids)
a <- fread("GSE92742_Broad_LINCS_sig_info.txt")
b <- filter(a, pert_type=='trt_sh.cgs'|pert_type=='trt_oe')
c <- cids[which(cids$sig_name %in% b$sig_id),]
index <- as.numeric(c$index)
d <- b[,c('sig_id','pert_iname')]
d1 <- merge(d,c,by.x='sig_id',by.y='sig_name',all.x=T)
d2 <- d1[order(d1$index),]

## map gene symbosl to ENSG
# ensembl = useEnsembl(biomart="ensembl", data="hsapiens_gene_ensembl")
# e <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart=ensembl)
# map <- merge(d1,e,by.x="pert_iname",by.y="hgnc_symbol",all.x=T)
# map_order <- map[order(map$index),]
# map_unique <- map_order[!duplicated(map_order$sig_id),]

##load data matrix of Z scores
ds <- parse.gctx("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", cid=index)

##map entrez genes to ENSG
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
f <- getBM(attributes=c('hgnc_symbol','entrezgene'), mart=ensembl)                
entrez <- as.data.frame(row.names(ds@mat))
colnames(entrez) <- c("entrez_id")
entrez$index <- row.names(entrez)
merge_entrez <- merge(entrez,f,by.x="entrez_id","entrezgene",all.x=T)
uniq_entrez <- merge_entrez[!duplicated(merge_entrez$entrez_id),]
sort_entrez <- uniq_entrez[order(as.numeric(uniq_entrez$index)),]

#colnames(ds@mat) <- map_unique$ensembl_gene_id

#distal <- as.data.frame(which(distant_mqtl$snp_hgnc_symbol %in% d2$pert_iname & distant_mqtl$cpg_hgnc_symbol %in% )
#distal1 <- distal[distal$s]

##only keep pairs which we can evaluate

##original
load("godmc_phase2_analysis/12_l1000-overlap/data/distant_mqtl.rdata")
distal <- as.data.table(distant_mqtl)
distal1 <- subset(distal, snp_hgnc_symbol %in% d2$pert_iname)
distal2 <- subset(distal1, cpg_hgnc_symbol %in% sort_entrez$hgnc_symbol)
  
datalist=list()
# pull out Z scores
for (i in 1:nrow(distal2)) {
  cols <- which(d2$pert_iname==as.character(distal2[i,3]))
  row <-which(sort_entrez$hgnc_symbol==as.character(distal2[i,4]))
  indices <- cbind(row,cols)
  #   if (nrow(d2[cols,])==1) {
  #     paste(distal2[i,],d2[cols,],ds@mat[indices])  
  #   } else {
  datalist[[i]] <- cbind(distal2[i,],d2[cols,][,1],ds@mat[indices])  
}
output = do.call(rbind, datalist)
colnames(output) <- c("cpg","snp","snp_hgnc_symbol","cpg_hgnc_symbol","sig_info","Z")
  
save(output, file = "perturbations_Zscores_perm.Rdata")

##permutations
for (j in 1:10) {
  load(paste("godmc_phase2_analysis/12_l1000-overlap/data/distant_perm", j, ".rdata", sep=""))
  distal <- as.data.table(perm)
  distal1 <- subset(distal, snp_hgnc_symbol %in% d2$pert_iname)
  distal2 <- subset(distal1, cpg_hgnc_symbol %in% sort_entrez$hgnc_symbol)

datalist=list()
# pull out Z scores
for (i in 1:nrow(distal2)) {
  cols <- which(d2$pert_iname==as.character(distal2[i,3]))
  row <-which(sort_entrez$hgnc_symbol==as.character(distal2[i,4]))
  indices <- cbind(row,cols)
#   if (nrow(d2[cols,])==1) {
#     paste(distal2[i,],d2[cols,],ds@mat[indices])  
#   } else {
  datalist[[i]] <- cbind(distal2[i,],d2[cols,][,1],ds@mat[indices])  
}
output = do.call(rbind, datalist)
colnames(output) <- c("cpg","snp","snp_hgnc_symbol","cpg_hgnc_symbol","sig_info","Z")

save(output, file = paste("perturbations_Zscores_perm", j, ".Rdata", sep=""))
}

# ds@mat[indices]
# indices$sig_id <- d2[cols,][,1]

##map sig names to ENSG
# a <- fread("GSE92742_Broad_LINCS_sig_info.txt")
# sig_names <- colnames(ds@mat)
# key <- a[which(sig_names %in% a$sig_id),]