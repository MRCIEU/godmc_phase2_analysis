#!/bin/bash

#PBS -N format-gwas-cc
#PBS -o /panfs/panasas01/shared-godmc/job_report/format-gwas-cc.o
#PBS -e /panfs/panasas01/shared-godmc/job_report/format-gwas-cc.e
#PBS -l walltime=20:00:00
#PBS -t 1-22
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash

set -e

echo "Running on ${HOSTNAME}"

chr=${PBS_ARRAYID}

# script to reformat GWAS p values

#The files must contain no header and have genomic
#position (build 37) in the first column and P-value from association analysis in the second column
#using a space as a column delimiter. 
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/09_garfield/scripts
JOBNO="`echo $PBS_JOBID | sed s/.bluequeue1.cvos.cluster//`" 
WORKDIR="/local/${PBS_O_LOGNAME}.${JOBNO}"
PVALDIR="/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/pval"
GWADIR="/projects/MRC-IEU/research/projects/ieu2/p5/021/working/data/Astle2016"
rsync -av epzjlm@bluecrystalp3.bris.ac.uk:$GWADIR/*gz $WORKDIR
cd $GWADIR
ls *-Build37.f.tsv.gz >bloodcounts.txt

filename="$1"
while read -r line;
do
name="$line"
name2=${name%-Build37.f.tsv.gz}
echo $name
echo $name2

mkdir -p $PVALDIR/$name2

#for chr in {1..22};
#      do 
      	zcat $name |awk -v i=$chr '$3==i{print $4, $10}'|sort -n >$PVALDIR/$name2/chr${chr};
#      	done

done <"bloodcounts.txt"


rm $WORKDIR