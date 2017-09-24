# Creating regulatory graphs using cis and trans mQTLs

Assume that cis-mQTLs can instrument the effects of one CpG on another. Find CpGs that have trans effects where the trans-SNP is a cis-mQTL for another CpG. 

```
Rscript creg_tcpg.r
```

This reads the `16_clumped.rdata` results to find all top trans hits. It then sequentially reads in the meta analysis results to identify which of the trans SNPs influence other CpGs in cis. These cis-trans putative associations are stored in `creg_tcpg.rdata`.

Next, use the coloc package to estimate the likelihood that the cis and trans mQTLs for a particular SNP share the same causal variant. 

```
sbatch coloc.sh
Rscript coloc_collate.r
```

This generates `coloc.rdata`. 

Merge the results from colocalisation analysis and construct a graph. Only keep SNPs that are Zhou filtered, have posterior probability of colocalisation > 0.8, trans pval < 1e-14 and cis pval < 1e-10.

Create a graph and use Walktrap algorithm to find communities of CpGs - these are CpGs that are connected to each other by some path in the graph.

```
Rscript make_graph.r
```

This generates `graph.rdata`.

Next things to do:

1. Is there enrichment of the genes within a community for GO terms etc?
2. If a CpG associates with a trait does that make other CpGs in the community more likely to associate with that trait?



Get annotations for 450k CpGs

```
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
```

