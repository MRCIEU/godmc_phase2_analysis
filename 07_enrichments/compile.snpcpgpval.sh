#!/bin/bash
  
#PBS -N compile   
#PBS -o compile-output
#PBS -e compile-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=16
#PBS -S /bin/bash
# PBS -t 1-10

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

dir="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16"
cd $dir
touch snpcpgpval.txt
touch tmp.txt
for i in `seq 1 962`;do
echo $i
zcat /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_$i.txt.gz |awk 'NR>1 {print $1,$8, '$i'}' >>snpcpgpval.txt
zcat /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_$i.txt.gz |awk -F'[ :_]' 'NR>1 {print $1,$2,$3}' >>tmp.txt
done

paste -d ' ' tmp.txt snpcpgpval.txt > snpcpgpval.tmp
mv snpcpgpval.tmp /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.txt
#rm /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/tmp.txt
gzip /panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/snpcpgpval.txt

for i in `seq 1 23`; do
echo $i
zcat $dir/snpcpgpval.txt.gz| awk -v mychr="chr$i" '$1==mychr {print $0}' > $dir/snpcpgpval.chr$i.txt
done
