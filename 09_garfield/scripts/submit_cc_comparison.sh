#!/bin/bash

#PBS -N cc
#PBS -o /panfs/panasas01/shared-godmc/job_report/cc.o
#PBS -e /panfs/panasas01/shared-godmc/job_report/cc.e
#PBS -l walltime=10:00:00
#PBS -t 1-36
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash
# PBS -q testq

set -e

echo "Running on ${HOSTNAME}"

i=${PBS_ARRAYID}
cd /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/09_garfield/scripts
echo $PBS_JOBID
JOBNO="`echo $PBS_JOBID | sed s/.master.cm.cluster//`" 
echo $JOBNO 
WORKDIR="/local/${PBS_O_LOGNAME}.${JOBNO}"
echo $WORKDIR
GWADIR="/projects/MRC-IEU/research/projects/ieu2/p5/021/working/data/Astle2016"
mkdir $WORKDIR


trait=`awk 'NR=='$i' {print $1}' <bloodcounts.txt`
echo $trait
traitname=${trait%-Build37.f.tsv.gz}
echo $traitname

rsync -av epzjlm@bluecrystalp3.bris.ac.uk:$GWADIR/$trait $WORKDIR
rsync -av epzjlm@bluecrystalp3.bris.ac.uk:$GWADIR/meta_data_forMRBase_Astle2016.txt $WORKDIR

cd $WORKDIR

zcat $trait |awk '{print $2,$3,$4,$5,$6,$9,$10,$11,$13}' > $traitname

Rscript /panfs/panasas01/sscm/epzjlm/repo/godmc_phase2_analysis/09_garfield/scripts/cc_comparison.R $traitname

rm -r $WORKDIR
