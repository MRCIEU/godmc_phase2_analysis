library(MSigDB)
library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other



# Download the list of gene names from http://www.genenames.org/cgi-bin/download
ok <- fread("zcat ../data/ok8ktmr.gz")

# Get the mappings of CpGs to gene names from IlluminaHumanMethylation450kanno.ilmn12.hg19
# Make sure that all names are current - check with previous names

masterlist <- data.frame(sym=ok[["Approved Symbol"]], prev=ok[["Previous Symbols"]], syn=ok$Synonyms, stringsAsFactors = FALSE)
masterlist$alt <- masterlist$prev
masterlist$alt[masterlist$prev == "" & masterlist$syn != ""] <- masterlist$syn[masterlist$prev == "" & masterlist$syn != ""]
masterlist$alt[masterlist$prev != "" & masterlist$syn != ""] <- paste(
	masterlist$prev[masterlist$prev != "" & masterlist$syn != ""],
	masterlist$syn[masterlist$prev != "" & masterlist$syn != ""], sep=", "
)
masterlist$alt <- gsub(", ", ",", masterlist$alt)

lst <- setNames(strsplit(as.character(masterlist$alt), ','), masterlist$sym)
m2 <- stack(lst)
m3 <- data.frame(values="", ind=masterlist$sym[masterlist$prev == "" & masterlist$syn == ""], stringsAsFactors = FALSE)
masterlist <- rbind(m2, m3)

anno2 <- data.frame(cpg=rownames(anno), sym=anno$UCSC_RefGene_Name, stringsAsFactors = FALSE)
anno3 <- stack(setNames(strsplit(anno2$sym, ";"), anno2$cpg))
anno3$ind <- as.character(anno3$ind)
anno4 <- data.frame(values="", ind=anno2$cpg[anno2$sym==""], stringsAsFactors = FALSE)
anno <- rbind(anno3, anno4)

annocur <- subset(anno, values %in% masterlist$ind)
annoold <- subset(anno, values %in% masterlist$values & values != "")
annoempty <- subset(anno, values == "")
ind <- match(annoold$values, masterlist$values)
annoold$values <- masterlist$ind[ind]

anno <- rbind(annocur, annoold, annoempty)


replace_old <- function(gl, masterlist, annolist)
{
	gl <- gl[(gl %in% masterlist$ind | gl %in% masterlist$values)]
	if(length(gl) > 0)
	{
		index <- gl %in% masterlist$values
		i <- match(gl[index], masterlist$values)
		gl[index] <- as.character(masterlist$ind[i])
	} else {
		message("none left")
		return(gl)
	}

	gl <- gl[gl %in% annolist]
	return(gl)
}

k <- 1
MSigDB2 <- list()
annolist <- unique(anno$values)
for(i in 1:length(MSigDB))
{
	MSigDB2[[names(MSigDB)[i]]] <- list()
	for(j in 1:length(MSigDB[[i]]))
	{
		message(k)
		out <- replace_old(MSigDB[[i]][[j]], masterlist, annolist)
		if(length(out) > 2)
		{
			MSigDB2[[i]][[names(MSigDB[[i]])[j]]] <- out
		}
		k <- k+1
	}
}


save(anno, masterlist, MSigDB2, file="../data/entrez_genes.rdata")

