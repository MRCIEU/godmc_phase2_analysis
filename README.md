# Phase 2 analysis

RDSF location:

```
/projects/MRC-IEU/research/data/godmc/_devs/GODMC_Analysis/data
```

BC4 location:

```
/mnt/storage/private/mrcieu/research/GODMC_Analysis/
```

BC3 location:

```
/panfs/panasas01/shared-godmc/
```

SFTP server location:

```
/srv/sftponly/GoDMC/
```

---

The paper is being written [here](https://drive.google.com/drive/folders/0B0vAR1k68I5fZkhpS1I3and0T2s?usp=sharing)

**Current strategy**

1. Use Google docs so that we can write and edit collaboratively
2. Use [paperpile](http://paperpile.com) for citations. Install it on google chrome so that you can edit / add citations

Problems with google docs: It doesn't do referencing of figures/tables automatically, so have to manually type e.g. Supplementary figure 1 etc.

The paper structure is currently written in the Main document in the shared google drive folder

## To do

There is a task list [here](https://docs.google.com/spreadsheets/d/1VihsoQhCNYwY07g-Asjgr6p9PBnQ44PQqhCSgzYm5mQ/edit?usp=sharing).

Information required from cohorts:

1. Informed consent information
2. Ethical approval information
3. Authorships


QC section for supplementary information:

1. Generate table from cohort_descriptives.rdata (03/analysis.r)
2. Plot of n CpGs and n SNPs (03/analysis.r)
3. M statistics (01/mstat.R)
4. comparison of FE and MRE (03/analysis.r)
5. Histograms of effect directions in meta analysis (03/analysis.r)
6. Meta regression (01/mstat.metaregression.R)
7. Check 16a and lambdas - add lambdas to table
8. Show that largest effect sizes have largest heterogeneity because standard errors are very small - scale issue


Main figures:

???

---

## How to run each section

## 01

This performs the meta analysis of the putative mQTLs from section 16 of the pipeline.

Each `cohort_16.tgz` file contains 962 `results.gz` files, which are in GWAMA format. All associations in each sub file are the same across cohorts, so to meta analyse we just need to extract the files to the right location, create the `metal` script, and run.

It is currently setup such that the `cohort_16.tgz` files need to be put into a directory somewhere and then gunzipped. Note - not untarred, e.g.

```
gunzip cohort_16.tgz
```

to produce `cohort_16.tar`. The reason for this is that the script extracts files from that tar file as it is needed, which is faster than doing it from a tgz file, and more space efficient than extracting everything.

Once you have a directory with all the `cohort_16.tar` files, then run

```
cd 01_meta_analysis_16
./run_metal.sh <batch number>
```

or 

```
qsub run_metal.sh
```

Explanation of the output:

```
# Fixed effects
Effect
StdErr
pval

# Heterogeneity stats
Direction
HetISq
HetChiSq # This is Q statistic
HetDf
HetPVal
tausq

# Additive random effects (DerSimonian-Laird estimator)
EffectARE
StdErrARE
PvalueARE

# Multiplicative random effects
# The Effect size for this is the same as the Fixed effects
StdErrMRE
PvalueMRE
```



## 02

Same as in 01, except the files sent by the cohorts are much larger (e.g. 20-30gb each). Each sub result file is in a specific binary format that is very space efficient. So the script extracts those files and runs an R script to convert the binary file into GWAMA format. The process continues as previously described from here.

Again, need to `gunzip` the `cohort_17.tgz` files to create `cohort_17.tar` files, and then

```
cd 02_meta_analysis_17
./run_metal.sh <batch number>
```

or 

```
qsub run_metal.sh
```

## 03

Performing clumping on the results of 01. Using a p-val threshold of 1e-4 for cis and 5e-8 for trans. Using a radius of 1Mb from CpG to denote cis/trans. This uses the European samples in the 1000 genomes data as an LD reference panel.

To run:

```
cd 03_clumping_16
./run_clumping.sh <batch number>
```

or

```
qsub run_clumping.sh
```

Once it is finished run 

```
Rscript aggregate_clumps.r
```

to create the file `results/16/16_clumped.rdata` which contains all the clumped results

To run this on bluecrystal 4 use

```
sbatch run_clumped_bc4.sh
```


Note: Get annotations for 450k CpGs

```r
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
```



## 04

Performing conditional analysis on the results of 01. Same p-val thresholds as in 03. Using a larger reference sample (4000 samples in ALSPAC, see `make_alspac_reference_dataset.sh` on details of how to make it.)

To run on bc4

```
cd 04_conditional_16
./run_conditional_bc4.sh <batch number>
```

or

```
sbatch run_conditional_bc4.sh
```



