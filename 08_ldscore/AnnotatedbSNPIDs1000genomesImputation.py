### from 1000 genomes format map to rsID in dbSNP 141

data = file("/mnt/data1/Eilis/Projects/References/SNPs/UCSC_All_SNPs_141_MapSNPIDs_UniqueMapping", "r")
data = data.readlines()

snpIDs = {}
for line in data:
    line = line.strip().split()
    if line[4] == "single":
        alleles = ["A", "T", "C", "G"]
    elif line[4] == "insertion":
	alleles = ["I", "R"]
    elif line[4] == "deletion":
	alleles = ["D", "R"]
    else: 
	alleles = ["NA", "NA"]
    if line[0].split("chr")[1] + ":" + line[1] + "_" + alleles[0] in snpIDs:
        for each in alleles:
            tmp = snpIDs[line[0].split("chr")[1] + ":" + line[1] + "_" + each]
            snpIDs[line[0].split("chr")[1] + ":" + line[1] + "_" + each] = tmp + "|" + line[2]

    else:
        for each in alleles:
            snpIDs[line[0].split("chr")[1] + ":" + line[1] + "_" + each] = line[2]
#    print line[0].split("chr")[1] + ":" + line[1] + "_" + alleles[1]


final = file("/mnt/data1/goDMC_Phase2/17_1_dbSNP141.txt", "w")
snps = file("/mnt/data1/goDMC_Phase2/17_1.txt", "r")
snps = snps.readlines()
final.write(snps[0] + "\trsID\tCpG\n")

for line in snps[1:]:
   
  line = line.strip().split()
  print line
  entries = line[0].split(":")
  print entries
  print entries[0][3:] +":" + entries[1] + "_" + line[1].upper()
  if entries[0][3:] +":" + entries[1] + "_" + line[1].upper() in snpIDs:
 
     name = snpIDs[entries[0][3:] +":" + entries[1] + "_" + line[1].upper()]
     cpg = entries[2].split("_")[1]
  else:
     name = entries[0][3:] + ":" + entries[1]
     cpg = entries[2].split("_")[1]
  line.append(name)
  line.append(cpg)
  # print line
  final.write("\t".join(line) + "\n")
 

