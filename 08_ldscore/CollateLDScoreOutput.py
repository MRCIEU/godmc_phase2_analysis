import os

final = file("LDScore_MergedOutput.txt", "w")
final.write("Probe\tnSNPswithSUmStats\tnSNPsAfterMergeWithRefPanel\tnSNPswithRegressionSNPLD\th2\tse\tLambdaGC\tMeanChi^2\tIntercept\tse\tRatio\n")


filesToRead = os.listdir("LDScore")

for each in filesToRead:
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
	
    final.write("\t".join(log) + "\n")

	
