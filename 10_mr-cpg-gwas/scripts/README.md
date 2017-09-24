# MR of CpGs on traits

- If a CpG only has a single SNP then use colocalisation analysis
- If a CpG has multiple conditionally independent SNPs then use IVW and Egger

First get the GWAS data:

```
Rscript cp_gwas_data.r
```

This copies the data from RDSF and is pretty slow. The summary datasets selected are listed in `ids.txt`. The data is stored in `../../data/gwas` with a file there called `00info.csv` that describes what each dataset is.

Next extract the clumped SNPs and conditionally independent SNPs from the GWAS datasets:

```
Rscript extract_clumped_from_gwas.r
```

This stores intermediate files in `../data/extracted` ready for conditional and colocalisation analysis.


## Colocalisation

Colocalisation analysis is only run on mQTLs that have p < 1e-5 in a particular GWAS. To setup the data:

```
Rscript setup_coloc.r
```

which will generate the necessary files in the correct format. Then submit:

```
sbatch coloc.sh
```

Currently 191 chunks, but may change as new data comes in and more putative associations found.

## MR

For CpGs that have more than one conditionally independent SNP can run MR. Use MendelianRandomization package because it handles correlated SNPs (IVW and Egger only).

Setup the data required:

```
Rscript setup_mr_ld.r
```

Run the analysis in 701 jobs (one for each GWAS): 

```
sbatch mr_ld.sh
```


