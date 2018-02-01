#!/bin/bash

# Parts that will stay the same throughout all clinical traits
DATADIR=/mnt/data1/programs/garfield-data
ANNOTDIR=/mnt/data1/goDMC_Phase2/godmc_phase2_analysis/09_garfield/annotationCisAmbvTrans
GARFIELDDIR=/mnt/data1/programs/garfield-v2

PRUNETAGSDIR=$DATADIR/tags/r01
CLUMPTAGSDIR=$DATADIR/tags/r08
MAFTSSDDIR=$DATADIR/maftssd_gc

ANNOTLINKFILE=$ANNOTDIR/link_file.txt
#NPERM=100000
NANNOT=$((`cat $ANNOTLINKFILE | wc -l` - 1))
PTHRESH=5e-8
BINNING=m5,n5,t5,p5,k5
CONDITION=0

# Parts that will change for each clinical trait
GWASTRAIT=$1
INPUTNAME=$2
PVALDIR=$DATADIR/pval/$INPUTNAME

OUTDIR=/mnt/data1/goDMC_Phase2/godmc_phase2_analysis/09_garfield/resultsCisAmbvTrans/$INPUTNAME
mkdir -p $OUTDIR

F1=$OUTDIR/garfield.prep.$INPUTNAME.out
F0=$OUTDIR/garfield.Meff.$INPUTNAME.out

cd $GARFIELDDIR
echo 'Prune and Clump'
echo -n > $F1
for CHR in `seq 1 22` #X
do
	echo 'CHR'$CHR
	./garfield-prep-chr -ptags $PRUNETAGSDIR/chr$CHR -ctags $CLUMPTAGSDIR/chr$CHR -maftss $MAFTSSDDIR/chr$CHR -pval $PVALDIR/chr$CHR -ann $ANNOTDIR/chr$CHR  -chr $CHR -o $F1 || { echo 'Failure!'; } 
done

echo 'Calculate effective number of annotations'
Rscript garfield-Meff-Padj.R -i $F1 -o $F0
NEA=$(head -1 $F0 |awk '{print $2}')
Padj=$(tail -1 $F0 |awk '{print $2}')

echo 'Calculate Enrichment and Significance'
F2=$OUTDIR/garfield.test.$INPUTNAME.out
Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -c $CONDITION
echo 'GARFIELD single annotation analysis complete'

echo 'Create Plots'
#Rscript garfield-plot.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

echo 'GARFIELD Analysis Complete!'
