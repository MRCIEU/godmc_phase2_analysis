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


#cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_mqtl
#if [ -e mqtlinputfile.txt ]
#then
#rm mqtlinputfile.txt
#fi

#touch mqtlinputfile.txt

#for i in `seq 1 22`; do
#echo $i
#cat chr$i >> mqtlinputfile.txt
#done

#cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_transmqtl

#if [ -e mqtlinputfile.txt ]
#then
#rm mqtlinputfile.txt
#fi

#touch mqtlinputfile.txt

#for i in `seq 1 22`; do
#echo $i
#cat chr$i >> mqtlinputfile.txt
#done


cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_mqtl_ihs
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_transmqtl_ihs
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_mqtl_fst
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_transmqtl_fst
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_mqtl_xpehhchb
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_transmqtl_xpehhchb
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_mqtl_xpehhyri
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/prep_transmqtl_xpehhyri
if [ -e mqtlinputfile.txt ]
then
rm mqtlinputfile.txt
fi

touch mqtlinputfile.txt

for i in `seq 1 22`; do
echo $i
cat chr$i >> mqtlinputfile.txt
done

cd /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield

mv /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/annotation/link_file_selection.txt /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield-data/annotation/link_file.txt

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_ihs/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/mqtl_ihs -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_ihs/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/transmqtl_ihs -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_fst/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/mqtl_fst -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_fst/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/transmqtl_fst -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_xpehhchb/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/mqtl_xpehhchb -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_xpehhchb/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/transmqtl_xpehhchb -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_mqtl_xpehhyri/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/mqtl_xpehhyri -g -m -t 0.0001 -progress

garfield-perm -n 100000 -a 1 -p 1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14 -pt 1e-10,1e-11,1e-12,1e-13,1e-14 -q n7,m10,t10 -i ../garfield-data/prep_transmqtl_xpehhyri/mqtlinputfile.txt \
-o /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/transmqtl_xpehhyri -g -m -t 0.0001 -progress


#Rscript garfield-plot.R /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/mqtl_iHS 100000 /panfs/panasas01/sscm/epzjlm/GARFIELD/garfield/results/figures/mqtl mqtl 5 0.05