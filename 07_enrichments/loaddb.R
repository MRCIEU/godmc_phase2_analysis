library(data.table)
library(LOLA)

#specify where to create the new LOLA db (use a temp folder, sothat the incomplete db doesn't interfere with usage)
wd="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states/25states"
setwd(wd)

#initialize and check the db
db=loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states/",collection="25states")

mkdir /panfs/panasas01/shared-godmc/godmc_phase2_analysis/gene_annotation/7regions
