# Phase 2 analysis

Note that the meta analysis, clumped and conditional analysis results from section 16 have now been uploaded to the SFTP server. These are stored in `/srv/sftponly/GoDMC/shared/16.tar`. To incorporate into your repository:

1. Update repository with `git pull` and navigate to the `results/` folder
2. Delete or `mv 16/ 16_old`
2. SFTP into server `sftp username@filetrn-scmv-d0.epi.bris.ac.uk`
3. Download the results `get shared/16.tar`
4. Untar `tar xvf 16.tar`
5. If you type `git status` it should not show any staged changes. Main files in here are `16_clumped.rdata`, `16_conditional.rdata` and the meta analysis results `16_*.txt.gz`

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


See documentation within each of the subsequent sections for details.


