#!/bin/bash

id=$1

dir="/mnt/storage/private/mrcieu/research/mr-eve"

Rscript extract_gwas.r \
--entities ../data/entity_info.rdata \
--bcf-dir $dir/gwas-files \
--out ../data/extract_gwas/$id.rdata \
--bfile $dir/vcf-reference-datasets/ukb/ukb_ref \
--vcf-ref $dir/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf \
--gwas-id $id \
--get-proxies yes

