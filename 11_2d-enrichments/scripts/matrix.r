library(dplyr)
library(doParallel)

load("../data/annotations.rdata")
load("../data/trans_clumped.rdata")

jid <- as.numeric(commandArgs(T)[1])
out <- paste0("../results/matrix/m", jid, ".rdata")

set.seed(jid)
if(jid != 0)
{
	message("Permuting")
	clumped$cpg <- sample(clumped$cpg)
}

an <- intersect(unique(cpgres$anno), unique(snpres$anno))
cpgres <- subset(cpgres, anno %in% an)
snpres <- subset(snpres, anno %in% an)

mat <- matrix(0, nrow(anno), nrow(anno))

(no_cores <- detectCores() - 1)
registerDoParallel(cores=no_cores)
pcl <- makeCluster(no_cores, type="FORK")

mat <- parLapply(pcl, 1:length(an), function(i)
{
	snps <- unique(subset(snpres, anno == an[i])$snp)
	cl <- subset(clumped, snp %in% snps)
	cpgs <- unique(cl$cpg)
	temp <- subset(cpgres, cpg %in% cpgs)
	temp$id <- 1:nrow(temp)
	cl2 <- merge(cl, temp, by="cpg")
	count <- group_by(cl2, anno) %>% summarise(n=n())
	out <- rep(0, nrow(anno))
	out[count$anno] <- count$n
	return(out)
})
stopCluster(pcl)

mat <- do.call(rbind, mat)
mat <- mat[,an]

save(mat, file=out)
