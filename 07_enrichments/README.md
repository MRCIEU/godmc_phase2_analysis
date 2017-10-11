`ld_regions.r` generates 

../results/enrichments/snpcontrolsets.rdata

- snpchr
- min
- max
- nproxies
- snppos
- MAF
- snptype
- SNP
- chromosome
- position
- closest450kcpg
- closest450kdistance
- tssdist
- closesttss
- cismQTL
- transmQTL
- mQTL

addcg_cpgfrq.R will add 

- GCcontent
- CpGcontent

../results/enrichments/snpcontrolsetsGC_CpGcontent.rdata

`dat` object

Each row represents an mQTL - association between a SNP and a CpG

- `cpgchr` = CpG chromosome
- `cpgpos` = CpG position
- `cpgname` = CpG Name
- `snpchr` = SNP chr
- `snppos` = SNP pos
- `snpname` = SNP name
- `effect_allele` = Effect allele
- `other_allele` = Other allele
- `effect_allele_freq` = Effect allele frequency
- `mqtl_effect` = Effect size of mQTL
- `mqtl_pval` = p-value (fixed effects meta analysis)
- `meta_directions` = For each of the cohorts, what was the direction of the effect size. ? = cohort did not have data
- `Isq` = Heterogeneity measure for meta analysis
- `samplesize` = Sample size
- `cis` = If true then CpG and SNP within 1  Mb
- `ld80_start` = The SNP position furthest upstream with LD rsq >= 0.8
- `ld80_end` = The SNP position furthest downstream with LD rsq >= 0.8
- `ld80_proxies` = Number of SNPs with LD rsq >= 0.8
- `ld80_dist` = Region size of SNPs with LD rsq >= 0.8
- `pchic` = Not NA if this mQTL has a SNP PCHIC chromosome interaction based on annotations from http://dx.doi.org/10.1016/j.cell.2016.09.037. The values represent the list of blood cell-types for which there was a significant chromosome interaction (chicago score > 5)

`cpg_clusters` object:

This is a list of CpGs that were identified to be in causal networks, where Mendelian randomisation was used to estimate if one CpG causally relates to another CpG. The cluster column is a numeric ID representing cluster membership.
