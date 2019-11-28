library(meffil)
path="/projects/MRC-IEU/research/data/aries/_devs/ARIES_Methylation_Tissue/data/"

samples<-read.table("~/LUMCsamples/samplesheet.GEO.txt",sep="\t",he=T)
names(samples)<-gsub("sex","Sex",names(samples))
names(samples)<-gsub("sentrix_id","Sample_Name",names(samples))
slide<-do.call("rbind",strsplit(as.character(samples$Sample_Name),split="_"))
samples<-data.frame(samples,slide=slide[,1])

samples2<-read.table("~/LUMCsamples/CORRECTED.csv",sep=",",he=T)
samples2$location<-gsub("C","C0",samples2$location)
samples2$location<-gsub("R","R0",samples2$location)
samples2$Sample_Name<-paste0(samples2$sentrix.id,"_",samples2$location)

samples$Basename<-paste0(path,samples$slide,"/",samples$Sample_Name)
qc.objects <- meffil.qc(samples, cell.type.reference="blood gse35069 complete", verbose=TRUE)

save(qc.objects,file="qc.objects.Robj")
length(qc.objects)
#161

qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5,
	snp.concordance.threshold             = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(
	qc.objects,
	parameters = qc.parameters,
	genotypes=NULL
)

save(qc.summary, file="qcsummary.Robj")

meffil.qc.report(qc.summary, output.file="qc-report.html")

outlier <- qc.summary$bad.samples
table(outlier$issue)
index <- outlier$issue %in% c("Control probe (dye.bias)", 
                              "Methylated vs Unmethylated",
                              "X-Y ratio outlier",
                              "Low bead numbers",
                              "Detection p-value",
                              "Sex mismatch",
                              "Genotype mismatch",
                              "Control probe (bisulfite1)",
                              "Control probe (bisulfite2)")

outlier <- outlier[index,]


outlier
#                        sample.name                      issue
#9305651099_R03C01 9305651099_R03C01 Methylated vs Unmethylated


y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="pc.fit.pdf",height=6,width=6)

pc <- 15
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste("norm.obj.pc",pc,".Robj",sep=""))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta,file=paste("norm.beta.pc",pc,".Robj",sep=""))
