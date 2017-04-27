# Phase 2 analysis

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



