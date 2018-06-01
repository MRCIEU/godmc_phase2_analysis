args <- commandArgs(T)
i <- as.numeric(args[1])

library(dplyr)

extract_from_17 <- function(fn, togrep)
{
	require(data.table)
	require(tidyr)
	message("Extracting")
	tmp <- tempfile() %>% basename
	outlist <- paste0(tmp, ".snplist")
	fnn <- paste0(tmp)
	write.table(togrep, file=outlist, row=F, col=F, qu=F)
	cmd <- paste0("zfgrep -f ", outlist, " ", fn, " > ", fnn)
	system(cmd)
	message("Reading")
	b <- fread(fnn)
	unlink(fnn)
	unlink(outlist)
	nom <- as.character(unlist(read.table(fn, nrows=1)))
	names(b) <- nom
	b$Pvalue <- as.numeric(b$Pvalue)
	b <- separate(b, MarkerName, c("snp", "cpg"), "_")
	return(b)
}

load("../results/mrbase_sig.rdata")
load("../data/snps_gwas.rdata")

chunks <- unique(sig$chunk)
# for(i in 1:length(chunks))
# {
	x <- subset(sig, chunk == chunks[i])
	y <- group_by(x, id.exposure) %>%
		do({
			y <- .
			temp <- subset(a, id.exposure == y$id.exposure)
			temp2 <- merge(y, temp, by="id.exposure", all=TRUE)
			data.frame(x=paste0(temp2$id, "_", temp2$outcome), stringsAsFactors = FALSE)
		})
	fn <- paste0("../../results/17/17_", chunks[i], ".txt.gz")
	exp <- extract_from_17(fn, y$x)
	save(exp, file=paste0("../results/mrbase_dat/dat", chunks[i], ".rdata"))
# }


