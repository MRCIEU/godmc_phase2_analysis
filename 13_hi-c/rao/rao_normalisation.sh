#!/bin/bash

#PBS -N hi-c-Rao_clean
#PBS -o /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/norm_clean-output
#PBS -e /panfs/panasas01/sscm/epwkb/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/job_reports/norm_clean-error
#PBS -l walltime=200:00:00
#PBS -l nodes=1:ppn=2
#PBS -t 1-100
#PBS -S /bin/bash

#1-253 chr_list

## NOTE: Submit 50 jobs to each array to prevent no space on device error: 
## some chr[1, 2 & 3] + chr[i] interactions take >200 hours to normalise 

set -e

echo $PBS_ARRAYID

echo "Running on ${HOSTNAME}"
start_time=`date +%s`

# Data
#*.RAWobserved = raw observed contact matrix in sparse matrix notation (mxn, more zero values than non-zero values)
# i, j, M_ij (value) e.g 100kb first row, 10th col of matrix = i=0 & j=900000 (only upper half is provided)
# All analyses and results in main paper were performed with KR normalised contact matrices
# KR norm vectors = 1st line = 1st row/col of the raw contact matrix and so on
# To normalise M_ij divide the entry by the corresponding norm factors for i and j
# e.g RAWobserved file:
# 40000 41000 59
# divide 59 by ((40000/1000)+1=41) 41st line and ((41000/1000)+1=42) 42nd line
# therefore = 59/(line_i*line_j)
# KR vectors i=first chr (left) and j=second chr (right)
# All data 1kb matrix resolution and > MAPQE30


# 1. Generate file of line numbers for norm
# 2. Split i/j chr into two vectors i/j (line per contact matrix)
# 3. Gen file of norm values for i/j
# 4. Gen norm values from RAW data (col3/i*j)

# File structure

#dir=("chr1_chr10" "chr1_chr11" "chr1_chr12" "chr1_chr13" "chr1_chr14" "chr1_chr15" "chr1_chr16" "chr1_chr17" "chr1_chr18" "chr1_chr19" "chr1_chr2" "chr1_chr20" "chr1_chr21" "chr1_chr22" "chr1_chr3" "chr1_chr4" "chr1_chr5" "chr1_chr6" "chr1_chr7" "chr1_chr8" "chr1_chr9" "chr1_chrX" "chr10_chr11" "chr10_chr12" "chr10_chr13" "chr10_chr14" "chr10_chr15" "chr10_chr16" "chr10_chr17" "chr10_chr18" "chr10_chr19" "chr10_chr20" "chr10_chr21" "chr10_chr22" "chr10_chrX" "chr11_chr12" "chr11_chr13" "chr11_chr14" "chr11_chr15" "chr11_chr16" "chr11_chr17" "chr11_chr18" "chr11_chr19" "chr11_chr20" "chr11_chr21" "chr11_chr22" "chr11_chrX" "chr12_chr13" "chr12_chr14" "chr12_chr15" "chr12_chr16" "chr12_chr17" "chr12_chr18" "chr12_chr19" "chr12_chr20" "chr12_chr21" "chr12_chr22" "chr12_chrX" "chr13_chr14" "chr13_chr15" "chr13_chr16" "chr13_chr17" "chr13_chr18" "chr13_chr19" "chr13_chr20" "chr13_chr21" "chr13_chr22" "chr13_chrX" "chr14_chr15" "chr14_chr16" "chr14_chr17" "chr14_chr18" "chr14_chr19" "chr14_chr20" "chr14_chr21" "chr14_chr22" "chr14_chrX" "chr15_chr16" "chr15_chr17" "chr15_chr18" "chr15_chr19" "chr15_chr20" "chr15_chr21" "chr15_chr22" "chr15_chrX" "chr16_chr17" "chr16_chr18" "chr16_chr19" "chr16_chr20" "chr16_chr21" "chr16_chr22" "chr16_chrX" "chr17_chr18" "chr17_chr19" "chr17_chr20" "chr17_chr21" "chr17_chr22" "chr17_chrX" "chr18_chr19" "chr18_chr20" "chr18_chr21" "chr18_chr22" "chr18_chrX" "chr19_chr20" "chr19_chr21" "chr19_chr22" "chr19_chrX" "chr2_chr10" "chr2_chr11" "chr2_chr12" "chr2_chr13" "chr2_chr14" "chr2_chr15" "chr2_chr16" "chr2_chr17" "chr2_chr18" "chr2_chr19" "chr2_chr20" "chr2_chr21" "chr2_chr22" "chr2_chr3" "chr2_chr4" "chr2_chr5" "chr2_chr6" "chr2_chr7" "chr2_chr8" "chr2_chr9" "chr2_chrX" "chr20_chr21" "chr20_chr22" "chr20_chrX" "chr21_chr22" "chr21_chrX" "chr22_chrX" "chr3_chr10" "chr3_chr11" "chr3_chr12" "chr3_chr13" "chr3_chr14" "chr3_chr15" "chr3_chr16" "chr3_chr17" "chr3_chr18" "chr3_chr19" "chr3_chr20" "chr3_chr21" "chr3_chr22" "chr3_chr4" "chr3_chr5" "chr3_chr6" "chr3_chr7" "chr3_chr8" "chr3_chr9" "chr3_chrX" "chr4_chr10" "chr4_chr11" "chr4_chr12" "chr4_chr13" "chr4_chr14" "chr4_chr15" "chr4_chr16" "chr4_chr17" "chr4_chr18" "chr4_chr19" "chr4_chr20" "chr4_chr21" "chr4_chr22" "chr4_chr5" "chr4_chr6" "chr4_chr7" "chr4_chr8" "chr4_chr9" "chr4_chrX" "chr5_chr10" "chr5_chr11" "chr5_chr12" "chr5_chr13" "chr5_chr14" "chr5_chr15" "chr5_chr16" "chr5_chr17" "chr5_chr18" "chr5_chr19" "chr5_chr20" "chr5_chr21" "chr5_chr22" "chr5_chr6" "chr5_chr7" "chr5_chr8" "chr5_chr9" "chr5_chrX" "chr6_chr10" "chr6_chr11" "chr6_chr12" "chr6_chr13" "chr6_chr14" "chr6_chr15" "chr6_chr16" "chr6_chr17" "chr6_chr18" "chr6_chr19" "chr6_chr20" "chr6_chr21" "chr6_chr22" "chr6_chr7" "chr6_chr8" "chr6_chr9" "chr6_chrX" "chr7_chr10" "chr7_chr11" "chr7_chr12" "chr7_chr13" "chr7_chr14" "chr7_chr15" "chr7_chr16" "chr7_chr17" "chr7_chr18" "chr7_chr19" "chr7_chr20" "chr7_chr21" "chr7_chr22" "chr7_chr8" "chr7_chr9" "chr7_chrX" "chr8_chr10" "chr8_chr11" "chr8_chr12" "chr8_chr13" "chr8_chr14" "chr8_chr15" "chr8_chr16" "chr8_chr17" "chr8_chr18" "chr8_chr19" "chr8_chr20" "chr8_chr21" "chr8_chr22" "chr8_chr9" "chr8_chrX" "chr9_chr10" "chr9_chr11" "chr9_chr12" "chr9_chr13" "chr9_chr14" "chr9_chr15" "chr9_chr16" "chr9_chr17" "chr9_chr18" "chr9_chr19" "chr9_chr20" "chr9_chr21" "chr9_chr22" "chr9_chrX")

#dir=("chr1_chr2")

for NUM in $PBS_ARRAYID; do
dir=($(sed "${NUM}q;d" ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/chr_list.txt))
done

# 1-4
for x in ${dir[@]}; do
echo $x
cd ~/GoDMC_Analysis/Hi-C/Rao2014/GM12878_combined_interchromosomal/1kb_resolution_interchromosomal/$x/MAPQGE30/
pwd
rm -f norm.j

#sort -gk 1 *.RAWobserved > ${x}_1kb.RAWobservedi
#sort -gk 2 *.RAWobserved > ${x}_1kb.RAWobservedj

awk '{ print $1/1000+1 }' ${x}_1kb.RAWobservedi > lines_i
awk '{ print $2/1000+1 }' ${x}_1kb.RAWobservedj > lines_j

IFS=_ read var1 var2 <<< "$x"
echo $var1
echo $var2

for NUMi in $(cat lines_i); do
sed "${NUMi}q;d" ${var1}_1kb.KRnorm >>norm.i
done

rm -f lines_i

for NUMj in $(cat lines_j); do
sed "${NUMj}q;d" ${var2}_1kb.KRnorm >>norm.j
done

rm -f lines_j

wc -l norm.i norm.j *.RAWobservedi *.RAWobservedj
paste ${x}_1kb.RAWobservedi norm.i > tmpi
paste ${x}_1kb.RAWobservedj norm.j > tmpj

sort -gk 1 tmpi > tmpi2
sort -gk 1 tmpj > tmpj2

paste tmpi2 tmpj2 > tmp
awk '{ print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8 }' tmp > tmp2

# remove NaNs
grep "NaN" tmp2 | wc -l
sed -i '/NaN/d' tmp2
wc -l tmp2

awk '{ print $1, "\t", $2, "\t", $3/($4*$5)}' tmp2 >${x}_1kb.NORMobserved

wc -l norm.i norm.j *.RAWobserved* ${x}_1kb.NORMobserved tmp*
rm -f tmp* *.RAWobservedi *.RAWobservedj norm*

done

#########

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
