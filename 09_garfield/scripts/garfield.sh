#!/bin/bash

# Parts that will stay the same throughout all clinical traits
DATADIR=/mnt/data1/programs/garfield-v2/garfield-data
ANNODIR=/mnt/data1/goDMC_Phase2/godmc_phase2_analysis/09_garfield/annotation/
GARFIELDDIR=/mnt/data1/programs/garfield-v2/

PRUNETAGSDIR=$DATADIR/tags/r01
CLUMPTAGSDIR=$DATADIR/tags/r08
MAFTSSDDIR=$DATADIR/maftssd
ANNOTDIR=$ANNODIR

ANNOTLINKFILE=$ANNOTDIR/link_file.txt
NPERM=100000
NANNOT=$((`cat $ANNOTLINKFILE | wc -l` - 1))
PTHRESH=5e-5,5e-6,5e-7,5e-8
PTHRESHTEST=5e-5,5e-6,5e-7,5e-8
BINNING=m10n10t7

# Parts that will change for each clinical trait
GWASTRAIT=$1
INPUTNAME=$2
PVALDIR=$DATADIR/pval/$INPUTNAME

OUTDIR=/mnt/data1/goDMC_Phase2/godmc_phase2_analysis/09_garfield/results/$INPUTNAME
mkdir -p $OUTDIR

F1=$OUTDIR/garfield.prep.$INPUTNAME.out
F2=$OUTDIR/garfield.perm.$INPUTNAME.out

echo 'Starting Trait:'$INPUTNAME

echo 'Prune and Clump'

echo -n > $F1
for CHR in `seq 1 22` X
do

echo 'CHR'$CHR

cd $GARFIELDDIR
./garfield-prep $PRUNETAGSDIR/chr$CHR $CLUMPTAGSDIR/chr$CHR $MAFTSSDDIR/chr$CHR $PVALDIR/chr$CHR $ANNOTDIR/chr$CHR 895,975,976,977,978,979,980 >> $F1 || { echo 'Failure!'; } 

done

echo 'Calculate Fold Enrichment and Significance'

# greedy permutation step
./garfield-perm -n $NPERM -a $NANNOT -p $PTHRESH -pt $PTHRESHTEST -q $BINNING -i $F1 -o $F2 -g -m 100 -t 0.0001

# original permutation step
#./garfield-perm -n $NPERM -a $NANNOT -p $PTHRESH -pt $PTHRESHTEST -q $BINNING -i $F1 -o $F2

#echo 'Create Plots'

#Rscript garfield-plot.R $F2 $NPERM $OUTDIR/garfield.${INPUTNAME} ${GWASTRAIT} 10 0

echo 'GARFIELD Analysis Complete!'
