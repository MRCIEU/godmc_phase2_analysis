#!/bin/bash
  
#PBS -N extract   
#PBS -o extract-output
#PBS -e extract-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash
# PBS -t 1-22

echo Running on host `hostname`
echo Time is `date`

mydir="/panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/12_external_proposals/CHARGE"
cd $mydir

#zcat ../../results/16/16_1.txt.gz |head -n1 >GlycaemicTraits_SNPs_38cohorts.txt
#for i in `seq 1 962`; do
#cat $i
#zcat ../../results/16/16_${i}.txt.gz |fgrep -f SNPs_GlycaemicTraits.formatted.txt >GlycaemicTraits_SNPs_38cohorts$i.txt
#cat GlycaemicTraits_SNPs_38cohorts$i.txt >>GlycaemicTraits_SNPs_38cohorts.txt
#rm GlycaemicTraits_SNPs_38cohorts$i.txt
#done

#cd /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/17

#zcat 17_1.txt.gz |head -n1 > $mydir/GlycaemicTraits_SNPs_38cohorts17.txt
#for i in `seq 1 300`; do
#cat $i
#zcat 17_${i}.txt.gz |fgrep -f $mydir/SNPs_GlycaemicTraits.formatted.txt > $mydir/GlycaemicTraits_SNPs_38cohorts17_$i.txt
#cat $mydir/GlycaemicTraits_SNPs_38cohorts17_$i.txt >> $mydir/GlycaemicTraits_SNPs_38cohorts17.txt
#rm GlycaemicTraits_SNPs_38cohorts17_$i.txt
#done
#for i in `seq 1 23`;do echo $i; zgrep chr${i}: <GlycaemicTraits_SNPs_38cohorts17.txt.gz >chr$i;done

R CMD BATCH filterSNPresults.R filterSNPresults.Rout
