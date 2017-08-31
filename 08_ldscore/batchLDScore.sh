#awk -F, '{if ($4 > 0.5) print $1}' /mnt/data1/450K_reference/vanDongen_HeritabilityEstimates.txt > /mnt/data1/goDMC_Phase2/godmc_phase2_analysis/08_ldscore/vanDongen_h2_0.5_cpgs.txt#

batch_number=${1}

ldsc="/mnt/data1/Eilis/Projects/LDScore/ldsc/"
refData="/mnt/data1/Eilis/Projects/LDScore/"

#mkdir -p SepFiles
#mkdir -p Munged
#mkdir -p LDScore


while read p; 
do 
   echo "Extracting results for " ${p}
   grep $p ../../17_${batch_number}_dbSNP141.txt | awk '{print $18, $17, $6, $2,$3, $8}' > SepFiles/tmp_$p.txt; 
   sed '1i SNP\tN\tbeta\tA1\tA2\tP' SepFiles/tmp_$p.txt > SepFiles/sumstats_$p.txt;
   rm SepFiles/tmp_$p.txt;
   python ${ldsc}/munge_sumstats.py --sumstats SepFiles/sumstats_$p.txt --out Munged/sumstats_$p --merge-alleles ${refData}/w_hm3.snplist;
   python ${ldsc}/ldsc.py --h2 Munged/sumstats_$p.sumstats.gz --ref-ld-chr ${refData}/resources/ --w-ld-chr ${refData}/resources/ --out LDScore/${p}_h2;
done < resources/cpglist_${batch_number}.txt

