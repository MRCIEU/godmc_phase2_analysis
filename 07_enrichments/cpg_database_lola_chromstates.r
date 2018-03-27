library(data.table)
library(LOLA)

#Use this script to split bed files by segmentation annotation (e.g. roadmap segmentation)


#1. download segmentation bed files:
# wget -r -A"*_segments.bed" --level 1 'http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/'

# split region files by segment

#specify where to create the new LOLA db (use a temp folder, sothat the incomplete db doesn't interfere with usage)
wd="/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states/25states"
setwd(wd)

#get annotation (downloaded from data source and modified/reduced as needed)
annot=fread("Consolidated_EpigenomeIDs_summary_Table_red.csv")

#list the downloaded combined bed files that are to be split
infiles=list.files(path="./egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/",pattern="*.bed",full.names=TRUE)

#function to split bed files and create file annotation
process_bed=function(file.name,seg,dt){
  write.table(dt,paste0(file.name,"_",seg,".bed"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  file_size=list(id=paste0(file.name,"_",seg) ,filename=paste0(file.name,"_",seg,".bed"),EID=sub("_.*","",file.name),size=length(dt[[1]]))
  return(file_size)
}

all_sizes=data.table()

for (file in infiles){
  cat(file,"\n")
  bed=fread(file)
  file.name=sub("\\.bed","",sub(".*/","",file))
  print(file.name)

  sizes=bed[,process_bed(file.name=file.name,seg=V4[1],dt=list(V1,V2,V3,V4)),by=V4]
  all_sizes=rbindlist(list(all_sizes,sizes) )

}
state<-read.table("/panfs/panasas01/shared-godmc/GARFIELDv2/garfield-data/segmentations/output_annotations/states.txt",he=T)
state$STATE<-paste0("E",state$STATE)
o<-order(as.character(state$STATE))
state<-state[o,]

#create index file (annotate segmentation code)
index_file=merge(all_sizes,annot,by="EID")

state$STATE
#[1] "E1"  "E10" "E11" "E12" "E13" "E14" "E15" "E16" "E17" "E18" "E19" "E2" 
#[13] "E20" "E21" "E22" "E23" "E24" "E25" "E3"  "E4"  "E5"  "E6"  "E7"  "E8" 
#[25] "E9" 
levels(as.factor(index_file$seg_code))
#[1] "E1"  "E10" "E11" "E12" "E13" "E14" "E15" "E16" "E17" "E18" "E19" "E2" 
#[13] "E20" "E21" "E22" "E23" "E24" "E25" "E3"  "E4"  "E5"  "E6"  "E7"  "E8" 
#[25] "E9" 

setnames(index_file,"V4","seg_code")
levels(as.factor(index_file$seg_code))
index_file[,_explanationseg:=state$NO.[as.factor(seg_code)],]
setcolorder(index_file,c("id","filename","EID","seg_code","seg_explanation" ,"size","ORDER","GROUP","Epigenome Mnemonic","Standardized Epigenome name","Epigenome name (from EDACC Release 9 directory)","ANATOMY","TYPE","LAB","AGE","SEX","ETHNICITY","Single Donor or Composite","Comments"))

write.table(index_file,"../index.txt",sep="\t",row.names=FALSE,quote=FALSE)


#initialize and check the db
db=loadRegionDB("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/chrom_states/25states")

##
