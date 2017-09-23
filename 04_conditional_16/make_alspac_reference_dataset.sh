# make the reference dataset for conditional analysis

cd ../data/ref

head -n 4000 /panfs/panasas01/shared/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/hrc/released/2015-10-30/data/derived/unrelated_ids/children/inclusion_list.txt > keeplist.txt

plink \
--bfile /panfs/panasas01/shared/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/filtered/bestguess/maf0.01_info0.8/combined/data \
--keep keeplist.txt \
--make-bed --out out

bfile="out"

echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original
awk '{if (($5 == "A" || $5 == "T" || $5 == "C" || $5=="G") &&  ($6 == "A" || $6 == "T" || $6 == "C" || $6=="G")) print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;else print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;}' ${bfile}.bim.original > ${bfile}.bim

SNPfail1="SNPfail1"

cp ${bfile}.bim ${bfile}.bim.original2
touch ${SNPfail1}
touch ${bfile}.duplicates.txt

Rscript ~/sandpit/godmc/relateds/resources/genetics/harmonization.R \
	${bfile}.bim \
	${SNPfail1}



cp ${bfile}.bim ${bfile}.bim.original3
awk '{
if (++dup[$2] > 1) { 
	print $1, $2".duplicate."dup[$2], $3, $4, $5, $6 
} else {
	print $0 
}}' ${bfile}.bim.original3 > ${bfile}.bim


