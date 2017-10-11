# L1000 overlaps with trans mQTL

Suppose that a trans mQTL signifies the change of gene expression at the SNP leading to the change of gene expression at the CpG. We can test empirically by looking at perturbation experiments that see the global gene expression changes that arise due to a particular gene expression being perturbed.

The L1000 platform provides such a resource

[https://www.biorxiv.org/content/early/2017/05/10/136168](https://www.biorxiv.org/content/early/2017/05/10/136168)

Website:

[https://clue.io/](https://clue.io/)

GEO repository:

[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)

This document describes the GEO repository

[https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#](https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#)

We are interested in the following file:

```
GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
```

Description of the results

> Q: What are features (rows) in the data matrix? Posed differently, what is the gene space accessible by L1000?
A: In L1000 datasets, features are genes and the matrix values correspond to their raw, normalized, or differential expression values, depending on which level of data is being used. L1000 reports on 12,328 unique genes; 978 of these are the landmark genes, which are directly measured. The remaining 11,350 are computationally inferred. 9,196 of these 11,350 genes are inferred with high fidelity, and together with the 978 landmarks comprise the Best INFerred Genes (BING) feature space, containing 10,174 genes total. We term the entire space of 12,328 genes as All Inferred Genes (AIG; see figure below).
>
> The unique identifier for each row is the Entrez ID for the gene.
>
> Note that earlier releases used Affymetrix-based identifiers. That is no longer necessary as data and inference is benchmarked against RNA-Seq datasets. Hence we use NCBI gene entrez gene identifiers (id and symbol).
>
> Because we now map to gene symbols, the number of inferred features in the current matrices provided is 12,328 (unique genes) and not 22,268 (which used to be the count based on earlier Affymetrix probe sets)

The gctx format can be accessed using the cmapR R package

[https://github.com/cmap/cmapR](https://github.com/cmap/cmapR)


---

There is also another similar experiment using shRNA (small-hairpin). There are only ~50 CRISPR modified genes in the first dataset, but in the shRNA dataset there are ~5000 modified genes.

[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)

Focus on this first, using distant mQTLs (where the SNP and CpG are on different chromosomes).

## Annotate SNP-CpG pairs

1. Use biomaRt to get entrez gene IDs and their ranges
2. Use CpG position to find if it lands within a gene range
3. Use LD mapping for each SNP (positions that span SNPs with LD >= 0.8) to find range for each SNP and find overlaps with gene ranges
4. For each SNP-CpG pair enumerate all possible combinations of SNP genes and CpG genes

To run this:

```
Rscript annotate_positions.r
```

## Find L1000 perturbations that overlap with SNP-CpG genes

...


## Test for significance

Permute the SNP-CpG list, then re-annotate



