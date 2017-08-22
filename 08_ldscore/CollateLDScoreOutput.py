import os

## load van Dongen SNP heritability estimates for comparision

snpH2 = file("/mnt/data1/450K_reference/vanDongen_HeritabilityEstimates.txt", "r")
snpH2 = snpH2.readlines()

probeH2 = {}

for line in snpH2:
   line = line.strip().split(",")
   probeH2[line[0]] = line[1:]

print "CpG heritabilities loaded"

final = file("LDScore_MergedOutput.txt", "w")
final.write("Probe\tnSNPswithSumStats\tnSNPsAfterMergeWithRefPanel\tnSNPswithRegressionSNPLD\th2\tse\tLambdaGC\tMeanChi^2\tIntercept\tse\tRatio\trMZ\trDZ\th2_twinAE\tc2_twinACE\td2_twinADE\tPval_C_ACE\tPval_D_ADE\th2_total\th2_SNPs\tPval_age_Ve\tPval_age_Va\tPval_sex_Va\tPval_sex_Ve\n")


filesToRead = os.listdir("LDScore")

for each in filesToRead:
    #print each
    probe = each.split("_")[0]
    data = file("LDScore/" + each, "r")
    data = data.readlines()
    log = [data[11].split("/")[1].split("_")[0]]
    log = log + [data[16].split()[4]]
    log = log + [data[22].split()[6]]
    log = log + [data[23].split()[6]]
    log = log + data[26].split()[4:]
    log = log + [data[27].split()[2]]
    log = log + [data[28].split()[2]]
    log = log + data[29].split()[1:]
    log = log + [data[30].split("(")[0]]
    
    if probe in probeH2:
        log = log + probeH2[probe]
    else:
        log = log + ["NA", "NA", "NA", "NA", "NA", "NA","NA", "NA", "NA","NA", "NA", "NA","NA"]
    
    log = [x.strip("(") for x in log]
    log = [x.strip(")")  for x in log]
	
    final.write("\t".join(log) + "\n")

	
