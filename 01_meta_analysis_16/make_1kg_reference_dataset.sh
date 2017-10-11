# make the reference dataset for conditional analysis

cd ../data/ref

bfile="eur"
out="eur2"

echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original
awk '{if (($5 == "A" || $5 == "T" || $5 == "C" || $5=="G") &&  ($6 == "A" || $6 == "T" || $6 == "C" || $6=="G")) print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;else print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;}' ${bfile}.bim.original > ${out}.bim

SNPfail1="SNPfail1"

cp ${out}.bim ${out}.bim.original2
touch ${SNPfail1}
touch ${out}.duplicates.txt

Rscript ~/sandpit/godmc/relateds/resources/genetics/harmonization.R \
	${out}.bim \
	${SNPfail1}



cp ${out}.bim ${out}.bim.original3
awk '{
if (++dup[$2] > 1) { 
	print $1, $2".duplicate."dup[$2], $3, $4, $5, $6 
} else {
	print $0 
}}' ${out}.bim.original3 > ${out}.bim

cp ${bfile}.bed ${out}.bed
cp ${bfile}.fam ${out}.fam

grep "duplicate" ${out}.bim | cut -f 2 >> ${SNPfail1}

plink --bfile ${out} --exclude ${SNPfail1} --make-bed --out ${out}
