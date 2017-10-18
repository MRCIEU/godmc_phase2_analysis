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

This generates `graph.rdata` and `graph_unpruned.rdata`. The unpruned one has lots of creg SNPs that are close together because it's not using clumped data. This is overcome simply in graph.rdata by, for each tcpg, finding all SNPs on the same chromosome and filtering based on 1mb radius. Dramatically reduces the number.

It also generates `grinfo.rdata` which summarises the clusters and is used extensively further on.

It also generates the GRanges objects required for LOLA analysis


## Are the CpGs that share causal variants correlated?

We can use ARIES data to estimate the correlations between the CpGs that share causal variants.

Run permutations of CpG pairs to obtain a null. Correlate the correlation value against the Wald ratio value. Though, note that this doesn't necessarily expect a large correlation because we don't assume that CpG1 causes CpG2.

```
Rscript community_correlations.r
```



## Gene set enrichment analysis using GSEA and MSigDB

This uses 6 databases of annotations. For each community it tests for enrichment/overlaps. Then it permutes the community-CpG matchings 1000 times to get null distributions. 

```
sbatch graph_annotations.sh
```


## Enrichments using LOLA

To run this first run `sbatch lola_cpg_core.sh`. After this, the other 3 lola scripts can be run.

```
sbatch lola_cpg_ext.sh
sbatch lola_snp_core.sh
sbatch lola_snp_ext.sh
```




Next things to do:

1. If a CpG associates with a trait does that make other CpGs in the community more likely to associate with that trait?


Get annotations for 450k CpGs

```
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
```


## Quick summary after LOLA

CpGs in communities tend to be at motifs related to cohesin - which relates to DNA organisation in the nucleosome.
SNPs that influence CpG communities tend to be in regions to do with histone marks or RNA polymerase. 
Pioneer transcription factors have a role in unwinding inaccessible DNA to initiate TF binding

H3K9me3
H3K4me1
H3K36me3

Use LD min/max for SNP region /panfs/panasas01/shared-godmc/1kg_reference_ph3/snpcontrolsets_selection.rdata
For each CpG enrichment motif, what are the SNPs enriched for
Add extended database which includes roadmap
Look at differences across tissues for the same motif



