load("results/cpg_trait_coloc.rdata")
a <- subset(res, nsnp >= 10 & H4 > 0.8 & p < 1e-10)
#b <- TwoSampleMR::available_outcomes()
b<-read.table("results/available_outcomes.txt",sep="\t",he=T)

c <- merge(a, b, by.x="outcome", by.y="id",all.x=T)
length(unique(c$exposure))
#1176
table(c$category)

#     Disease Immune system   Metabolites   Risk factor 
#          485             0           755           707 
c_dis<-c[which(c$category=="Disease"),]
dim(c_dis)
df_dis<-data.frame(table(c_dis$subcategory))
df_dis<-df_dis[df_dis$Freq>0,]

df<-data.frame(table(c$trait))
df2<-df[df$Freq>0,]
dim(df2)
#220

#df2[df2$Freq>25,]
#                           Var1 Freq
#328              Crohns disease  107
#595                      Height  229
#660  Inflammatory bowel disease   77
#1017       Rheumatoid arthritis  104
#1020              Schizophrenia  101

cat<-data.frame(table(c$subcategory))
cat<-cat[cat$Freq>0,]
o<-order(cat$Freq)
cat<-cat[o,]

#gwas-cpgs
load("../06_mr-gwas-cpg/results/tophits_followup.rdata")
a <- subset(res, nsnp >= 10 & pval < 1e-10)
length(unique(a$exposure))
df<-data.frame(table(a$exposure))

g<-grep("Multiple sclerosis",res$exposure)
length(res[g,"outcome"])
#[1] 2472
length(unique(res[g,"outcome"]))
#[1] 789

df<-data.frame(table(a$exposure))
#> df
#                                                           Var1  Freq
#1                 Acute lymphoblastic leukemia (childhood) (NA)   278  
#2            Adiponectin levels (unit increase) (unit increase)    10
#3                  Age at menarche || ReproGen || 2014 || years    29
#4                         Age-related macular degeneration (NA)     7
#5               Body mass index (unit increase) (unit increase)     1
#6                              Breast cancer (early onset) (NA)  2989
#7                                           Celiac disease (NA)   818
#8                      Celiac disease || NA || 2011 || log odds    36
#9            Cholesterol, total (unit increase) (unit increase)     5
#10            HDL cholesterol (mg/dL increase) (mg/dL increase)    42
#11              HDL cholesterol (s.d. increase) (s.d. increase)     1
#12                       Height (unit increase) (unit increase)     2
#13                       Idiopathic membranous nephropathy (NA) 38269
#14                                         IgA nephropathy (NA)     3
#15              Multiple sclerosis || IMSGC || 2011 || log odds    56
#16              Multiple sclerosis || IMSGC || 2013 || log odds   204
#17                                      Multiple sclerosis (NA)    85
#18                                     Parkinson's disease (NA)     8
#19              Phospholipid levels (plasma) (DPA) (% increase)   210
#20                               Rheumatoid arthritis (EA) (NA)   730
#21                                    Rheumatoid arthritis (NA)  2864
#22               Rheumatoid arthritis || NA || 2014 || log odds     1
#23                                           Schizophrenia (NA)     3
#24  Serum protein levels (sST2) (unit increase) (unit increase)    10
#25 Trans fatty acid levels (Cis/trans-18:2, EA) (unit increase) 87226
#26                  Triglycerides || GLGC || 2013 || SD (mg/dL)     3
#27                                         Type 1 diabetes (NA)    68
#28                                      Ulcerative colitis (NA)    18
#29                               Urate || GUGC || 2013 || mg/dl     1
#30               Urate levels (mg/dl increase) (mg/dl increase)    22

a$exposure_recoded<-a$exposure
w<-which(a$exposure%in%c("Multiple sclerosis || IMSGC || 2011 || log odds","Multiple sclerosis || IMSGC || 2013 || log odds","Multiple sclerosis (NA)"))
a$exposure_recoded[w]<-"Multiple sclerosis (NA)"
w<-which(a$exposure%in%c("Rheumatoid arthritis (EA) (NA)","Rheumatoid arthritis || NA || 2014 || log odds","Rheumatoid arthritis (NA)"))
a$exposure_recoded[w]<-"Rheumatoid arthritis (NA)"
w<-which(a$exposure%in%c("HDL cholesterol (mg/dL increase) (mg/dL increase)","HDL cholesterol (s.d. increase) (s.d. increase)"))
a$exposure_recoded[w]<-"HDL cholesterol (NA)"
w<-which(a$exposure%in%c("Urate || GUGC || 2013 || mg/dl","Urate levels (mg/dl increase) (mg/dl increase)"))
a$exposure_recoded[w]<-"Urate levels (NA)"
w<-which(a$exposure%in%c("Celiac disease (NA)","Celiac disease || NA || 2011 || log odds"))
a$exposure_recoded[w]<-"Celiac disease (NA)"
a$Trait<-gsub("\\(.*?\\)","",a$exposure_recoded)
a$Trait<-gsub("\\|.*","",a$Trait)
a$Trait<-sub("\\s+$", "", a$Trait)
a$Trait<-gsub("\\'", "", a$Trait)


m<-match(unique(a$Trait),as.character(b$trait))
b2<-data.frame(unique(a$Trait),b[m,c("trait","subcategory")])
w<-which(b2[,1]%in%c("Age-related macular degeneration"))
b2[w,"subcategory"]<-"Eye"
w<-which(b2[,1]%in%c("Cholesterol, total","Adiponectin levels","Trans fatty acid levels","Phospholipid levels"))
b2[w,"subcategory"]<-"Lipid"
w<-which(b2[,1]%in%c("Trans fatty acid levels"))
b2[w,"subcategory"]<-"Fatty acid"
w<-which(b2[,1]%in%c("Urate levels"))
b2[w,"subcategory"]<-"Other"

w<-which(b2[,1]%in%c("Type 1 diabetes"))
b2[w,"subcategory"]<-"Diabetes"
w<-which(b2[,1]%in%c("Breast cancer","Acute lymphoblastic leukemia"))
b2[w,"subcategory"]<-"Cancer"
w<-which(b2[,1]%in%c("Idiopathic membranous nephropathy"))
b2[w,"subcategory"]<-"Kidney"
w<-which(b2[,1]%in%c("Serum protein levels"))
b2[w,"subcategory"]<-"Protein"

m<-match(a$Trait,as.character(b2[,1]))
dfcat<-data.frame(table(b2[m,c("subcategory")]))
dfcat<-dfcat[dfcat$Freq>0,]

df3<-unique(data.frame(a$exposure_recoded,a$outcome))
df3<-data.frame(table(df3[,1]))

data.frame(table(a$exposure_recoded))
#                                                           Var1  Freq
#1                 Acute lymphoblastic leukemia (childhood) (NA)   278
#2            Adiponectin levels (unit increase) (unit increase)    10
#3                  Age at menarche || ReproGen || 2014 || years    29
#4                         Age-related macular degeneration (NA)     7
#5               Body mass index (unit increase) (unit increase)     1
#6                              Breast cancer (early onset) (NA)  2989
#7                                           Celiac disease (NA)   854
#8            Cholesterol, total (unit increase) (unit increase)     5
#9                                          HDL cholesterol (NA)    43
#10                       Height (unit increase) (unit increase)     2
#11                       Idiopathic membranous nephropathy (NA) 38269
#12                                         IgA nephropathy (NA)     3
#13                                      Multiple sclerosis (NA)   345
#14                                     Parkinson's disease (NA)     8
#15              Phospholipid levels (plasma) (DPA) (% increase)   210
#16                                    Rheumatoid arthritis (NA)  3595
#17                                           Schizophrenia (NA)     3
#18  Serum protein levels (sST2) (unit increase) (unit increase)    10
#19 Trans fatty acid levels (Cis/trans-18:2, EA) (unit increase) 87226
#20                  Triglycerides || GLGC || 2013 || SD (mg/dL)     3
#21                                         Type 1 diabetes (NA)    68
#22                                      Ulcerative colitis (NA)    18
#23                                            Urate levels (NA)    23


