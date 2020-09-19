# Simulations to evaluate impacts of sparsity on colocalisation

In the GoDMC analysis, the two stage mQTL detection process means that the entire region of a significant signal is not necessarily retained. It will favour SNPs with lower p-values, which are more likely to be those in LD with the causal variant.

Objective: Investigate if the performance of coloc is adversely affected by the sparseness of the GWAS summary data.

## Approach

Need to simulate GWAS summary data for a region where there is a GWAS signal, and where the marginal effect estimates of each recorded variant is based on LD with the causal variants in the region and some sampling variation.

The method for simulating effect sizes and standard errors with appropriate sampling error for a given sample size is:

https://explodecomputer.github.io/simulateGP/articles/gwas_summary_data.html

Most of the theory for simulating under an LD matrix is laid out in:

https://github.com/explodecomputer/simulateGP/blob/master/vignettes/gwas_summary_data_ld.Rmd

Use the ALSPAC LD reference panel ~8000 samples as has been used for other parts of the GoDMC analysis.

To determine the SNPs to retain in a region, sample some proportion with the chance of being retained being higher if the p-value is higher. This mimics what is done in GoDMC, where an imperfect predictor of a well performing SNP is selected in stage 1 and then tested in stage 2.

To deal with coverage being less than 100%, either only look at the available SNPs, or fill in the missing SNPs to have p-values of 1.

## Contrasts

Allow the following to vary

- Numbers of causal variants in the region for each trait
- Genomic regions
- Variance explained by causal variants in each trait
- Coverage of the region
- Methods for dealing with incomplete coverage
- Whether the traits colocalise at the region or not


## To run

```
snakemake
```

