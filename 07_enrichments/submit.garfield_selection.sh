#!/bin/bash
  
#PBS -N garfield   
#PBS -o garfield-output
#PBS -e garfield-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
# PBS -t 1-22

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_mqtl_ihs
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_transmqtl_ihs
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_mqtl_fst
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_transmqtl_fst
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_mqtl_xpehhchb
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_transmqtl_xpehhchb
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_mqtl_xpehhyri
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_transmqtl_xpehhyri
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_mqtl_sds
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/prep_transmqtl_sds
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
#awk -F ' ' '{print "chr"'$i'":"$1,$0}' <chr$i >chr$i.tmp
#mv chr$i.tmp chr$i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield

cp /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/annotation/link_file_selection.txt /panfs/panasas01/shared-godmc/GARFIELD/garfield-data/annotation/link_file.txt

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_ihs/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_ihs -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_ihs/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_ihs -g -m -t 0.0001 -progress

perl -pe 's/mqtls/transmqtl/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_ihs >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_ihs.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_ihs.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_ihs

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_fst/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_fst -g -m -t 0.0001 -progress

perl -pe 's/ihs/fst/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_fst >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_fst.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_fst.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_fst

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_fst/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_fst -g -m -t 0.0001 -progress

perl -pe 's/ihs/fst/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_fst | perl -pe 's/mqtls/transmqtl/g' >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_fst.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_fst.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_fst    

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_xpehhchb/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhchb -g -m -t 0.0001 -progress

perl -pe 's/ihs/xpehhchb/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhchb >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhchb.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhchb.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhchb    

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_xpehhchb/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhchb -g -m -t 0.0001 -progress

perl -pe 's/ihs/xpehhchb/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhchb | perl -pe 's/mqtls/transmqtl/g' >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhchb.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhchb.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhchb    

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_xpehhyri/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhyri -g -m -t 0.0001 -progress

perl -pe 's/ihs/xpehhyri/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhyri >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhyri.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhyri.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_xpehhyri    

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_xpehhyri/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhyri -g -m -t 0.0001 -progress

perl -pe 's/ihs/xpehhyri/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhyri | perl -pe 's/mqtls/transmqtl/g' >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhyri.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhyri.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_xpehhyri    

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_sds/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_sds -g -m -t 0.0001 -progress

perl -pe 's/ihs/sds/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_sds >/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_sds.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_sds.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/mqtl_sds    

garfield-perm -n 100000 -a 5 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_sds/mqtlinputfile.txt \
-o /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_sds -g -m -t 0.0001 -progress

perl -pe 's/ihs/sds/g' </panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_sds | perl -pe 's/mqtls/transmqtl/g'>/panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_sds.tmp
mv /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_sds.tmp /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/transmqtl_sds    

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield/results

if [ -e selectionmetrics.txt ]
then
rm selectionmetrics.txt
fi

head -n 1 mqtl_ihs >selectionmetrics.txt
awk -F " " '$1==0 {print $0}' <mqtl_ihs |tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==0 {print $0}' <transmqtl_ihs | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==1 {print $0}' <mqtl_fst | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==1 {print $0}' <transmqtl_fst | tail -q -n +2  >>selectionmetrics.txt
awk -F " " '$1==2 {print $0}' <mqtl_xpehhchb | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==2 {print $0}' <transmqtl_xpehhchb | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==3 {print $0}' <mqtl_xpehhyri | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==3 {print $0}' <transmqtl_xpehhyri | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==4 {print $0}' <mqtl_sds | tail -q -n +2 >>selectionmetrics.txt
awk -F " " '$1==4 {print $0}' <transmqtl_sds | tail -q -n +2 >>selectionmetrics.txt

head -n1 selectionmetrics.txt >selection.txt
grep mqtls selectionmetrics.txt >>selection.txt 
head -n1 selectionmetrics.txt >selection_trans.txt
grep transmqtl selectionmetrics.txt >>selection_trans.txt 

cd /panfs/panasas01/shared-godmc/GARFIELD/garfield
Rscript garfield-plot.R /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/selection.txt 100000 /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/figures/mqtl mqtl 5 0.05
Rscript garfield-plot.R /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/selection_trans.txt 100000 /panfs/panasas01/shared-godmc/GARFIELD/garfield/results/figures/transmqtl transmqtl 5 0.05



