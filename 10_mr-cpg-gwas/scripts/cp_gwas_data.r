library(TwoSampleMR)
library(dplyr)
library(doParallel)


ids <- scan("ids.txt")
ao <- available_outcomes()
fn <- subset(ao, id %in% ids, select=c(path, filename, trait, sample_size, ncase, ncontrol, nsnp, id, category, sex))
fn <- fn[order(fn$sample_size, decreasing=TRUE), ]
fn <- subset(fn, !duplicated(trait))
fn$newfile <- paste0(gsub(" ", "_", paste(tolower(fn$category), tolower(gsub("[^[:alnum:] ]", "", fn$trait)), sep="__")), ".txt.gz")


###

dir.create("../../data/gwas", recursive=TRUE)
(no_cores <- detectCores() - 1)
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores, type="FORK")
result <- parLapply(cl, 1:nrow(fn), function(i)
{
	message(i)
	path <- paste0("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/passqc/", fn$filename[i])
	cmd <- paste0("gzip -c ", path, " > ../../data/gwas/", fn$newfile[i])
	system(cmd)
})
stopCluster(cl)
write.csv(fn, file="../../data/gwas/00info.csv")
