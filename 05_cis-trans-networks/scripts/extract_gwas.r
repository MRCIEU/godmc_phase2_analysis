if(!require(gwasvcftools))
{
	if(!required(devtools)) install.packages("devtools")
	devtools::install_github("MRCIEU/gwasvcftools")
}
library(gwasvcftools)
library(argparse)

# create parser object
parser <- ArgumentParser()
parser$add_argument('--entities', required=TRUE)
parser$add_argument('--bcf-dir', required=TRUE)
parser$add_argument('--gwas-id', required=TRUE)
parser$add_argument('--out', required=TRUE)
parser$add_argument('--bfile', required=TRUE)
parser$add_argument('--get-proxies', default='yes')
parser$add_argument('--vcf-ref', required=FALSE)
parser$add_argument('--tag-r2', type="double", default=0.6)
parser$add_argument('--tag-kb', type="double", default=5000)
parser$add_argument('--tag-nsnp', type="double", default=5000)
parser$add_argument('--palindrome-freq', type="double", default=0.4)
parser$add_argument('--threads', type="integer", default=1)
parser$add_argument('--no-clean', action="store_true", default=FALSE)

args <- parser$parse_args()
print(args)
tempname <- tempfile(pattern="extract", tmpdir=dirname(args[['out']]))
bcf <- file.path(args[['bcf_dir']], args[['gwas_id']], "harmonised.bcf")
load(args[['entities']])

snplist <- subset(entities, type != "creg_cpg") %$%
	data_frame(V1=snp_chr, V2=snp_pos, V3=snp_a1, V4=snp_a2, V5="b37", V6=snp_rsid) %>%
	subset(!duplicated(paste(V1,V2)))

extracted <- gwasvcftools::extract(bcf, snplist, tempname, args[['get_proxies']], args[["bfile"]], args[["vcf_ref"]], threads=args[["threads"]])
extracted$mrbaseid <- args[['gwas_id']]

save(extracted, file=args[["out"]])

