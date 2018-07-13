#echo 'Create Plots'

cd /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2

ANNOTNAME="ihs"
INPUTNAME=mqtl_${ANNOTNAME}

DIR="../garfield-data"
DATADIR="/panfs/panasas01/shared-godmc/GARFIELD/garfield-data"
ANNOTDIR=$DIR/annotation_selection
ANNOTLINKFILE=$ANNOTDIR/link_file_selection.txt
OUTDIR=$DIR/output

#Index Annotation Celltype Tissue Type Category
#0 ihs blood blood mqtls mqtl
#1 fst blood blood mqtls mqtl
#2 xpehhchb blood blood mqtls mqtl
#3 xpehhyri blood blood mqtls mqtl
#4 sds blood blood mqtls mqtl


cat $OUTDIR/mqtl_sds/header.txt $OUTDIR/mqtl_ihs/garfield.test.mqtl_ihs.out $OUTDIR/mqtl_fst/garfield.test.mqtl_fst.out  $OUTDIR/mqtl_xpehhchb/garfield.test.mqtl_xpehhchb.out $OUTDIR/mqtl_xpehhyri/garfield.test.mqtl_xpehhyri.out $OUTDIR/mqtl_sds/garfield.test.mqtl_sds.out>$OUTDIR/selection.txt
F2="$OUTDIR/selection.txt"

F0=$OUTDIR/$INPUTNAME/garfield.Meff.$INPUTNAME.out

Padj=$(tail -1 $F0 |awk '{print $2}')

Rscript /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2/garfield-plot_selection.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

#
ANNOTLINKFILE=$ANNOTDIR/link_file_selection_trans.txt
cat $OUTDIR/mqtl_sds/header.txt $OUTDIR/transmqtl_ihs/garfield.test.transmqtl_ihs.out $OUTDIR/transmqtl_fst/garfield.test.transmqtl_fst.out  $OUTDIR/transmqtl_xpehhchb/garfield.test.transmqtl_xpehhchb.out $OUTDIR/transmqtl_xpehhyri/garfield.test.transmqtl_xpehhyri.out $OUTDIR/transmqtl_sds/garfield.test.transmqtl_sds.out>$OUTDIR/selection_trans.txt
F2="$OUTDIR/selection_trans.txt"

F0=$OUTDIR/$INPUTNAME/garfield.Meff.$INPUTNAME.out
Padj=$(tail -1 $F0 |awk '{print $2}')

Rscript /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2/garfield-plot_selection.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

#
ANNOTLINKFILE=$ANNOTDIR/link_file_selection_cis.txt
cat $OUTDIR/mqtl_sds/header.txt $OUTDIR/cismqtl_ihs/garfield.test.cismqtl_ihs.out $OUTDIR/cismqtl_fst/garfield.test.cismqtl_fst.out  $OUTDIR/cismqtl_xpehhchb/garfield.test.cismqtl_xpehhchb.out $OUTDIR/cismqtl_xpehhyri/garfield.test.cismqtl_xpehhyri.out $OUTDIR/cismqtl_sds/garfield.test.cismqtl_sds.out>$OUTDIR/selection_cis.txt
F2="$OUTDIR/selection_cis.txt"

F0=$OUTDIR/$INPUTNAME/garfield.Meff.$INPUTNAME.out
Padj=$(tail -1 $F0 |awk '{print $2}')

Rscript /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2/garfield-plot_selection.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

#
ANNOTLINKFILE=$ANNOTDIR/link_file_selection_ambivalent.txt
cat $OUTDIR/mqtl_sds/header.txt $OUTDIR/ambivalentmqtl_ihs/garfield.test.ambivalentmqtl_ihs.out $OUTDIR/ambivalentmqtl_fst/garfield.test.ambivalentmqtl_fst.out  $OUTDIR/ambivalentmqtl_xpehhchb/garfield.test.ambivalentmqtl_xpehhchb.out $OUTDIR/ambivalentmqtl_xpehhyri/garfield.test.ambivalentmqtl_xpehhyri.out $OUTDIR/ambivalentmqtl_sds/garfield.test.ambivalentmqtl_sds.out>$OUTDIR/selection_ambivalent.txt
F2="$OUTDIR/selection_ambivalent.txt"

F0=$OUTDIR/$INPUTNAME/garfield.Meff.$INPUTNAME.out
Padj=$(tail -1 $F0 |awk '{print $2}')

Rscript /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2/garfield-plot_selection.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj


###

#perl -pe 's/mqtls/transmqtls/g' < $F2 >$F2.tmp
#cat $OUTDIR/selection.txt $F2.tmp > $OUTDIR/selection_all.txt


#ANNOTLINKFILE=$ANNOTDIR/link_file_selection_all.txt
#F2="$OUTDIR/selection_all.txt"
#Rscript /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2/garfield-plot_selection.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj


# echo 'Prioritize relevant annotations by conditional analysis'
# Additional prioritization of annotations
# CONDITION=1
# CONDITIONTHRESH=0.05
# Rscript garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -s $SUBSET -c $CONDITION -ct $CONDITIONTHRESH -padj $Padj
#	echo 'GARFIELD model selection complete'
