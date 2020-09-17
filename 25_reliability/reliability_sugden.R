library(data.table)
library(ggplot2)
library(dplyr)

retaincpg <- scan("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/retain_from_zhou.txt", what="character")
 
#exclusion probes from TwinsUK
excl<-read.table("~/repo/godmc_phase1_analysis/07.snp_cpg_selection/data/450k_exclusion_probes.txt",he=T)
#42446
rm<-which(retaincpg%in%excl[,1])
#14882
retaincpg<-retaincpg[-rm]

r<-fread("Sugden_MethylationReliability_Data_S1.txt",sep="\t",he=T)
r<-as.data.frame(r)
poor<-r[which(r$Reliability<0.4),]
w<-which(poor[,1]%in%retaincpg)
length(w) #299287
nrow(poor) #337550

good<-r[which(r$Reliability>0.75),]
nrow(good) #28249
i<-intersect(good[,1],retaincpg)
length(i) #19199

w<-which(r[,1]%in%retaincpg)
mean(r[w,"Reliability"])
#0.1961723
median(r[w,"Reliability"])

mean(r[-w,"Reliability"])
#0.3116608
median(r[-w,"Reliability"])

#w<-which(r[,1]%in%retaincpg)
#r$MASK_general<-"TRUE"
#r$MASK_general[w]<-"FALSE" #retain

mf <- fread("HM450.hg19.manifest.tsv.gz", header=TRUE)
mf<-as.data.frame(mf)
a<-mf[which(mf$MASK_general==FALSE),] #retain

mf2<-merge(r, mf[,c(5,47:57)],by.x=1,by.y="probeID")

load("/panfs/panasas01/shared-godmc/godmc_phase2_analysis/results/16/16_clumped.rdata")
clumped<-as.data.frame(clumped)
clumped2<-clumped[which(clumped$cis==TRUE & clumped$pval<1e-8 | clumped$cis==FALSE & clumped$pval<1e-14),]
mqtlcpg<-unique(clumped2$cpg)

w<-which(mf2[,1]%in%mqtlcpg)
mf2$mqtlprobe<-"FALSE"
mf2$mqtlprobe[w]<-"TRUE"

rt<-which(mf2[,1]%in%retaincpg)

mf2[rt,] %>%                                        # Specify data frame
  group_by(mqtlprobe) %>%                         # Specify group indicator
  summarise_at(vars(Reliability),              # Specify column
               list(name = median))               # Specify function

mf2 %>%                                        # Specify data frame
  group_by(MASK_typeINextBaseSwitch) %>%                         # Specify group indicator
  summarise_at(vars(Reliability),              # Specify column
               list(name = median))

 mf2 %>%                                        # Specify data frame
  group_by(MASK_snp5_GMAF1p) %>%                         # Specify group indicator
  summarise_at(vars(Reliability),              # Specify column
               list(name = median))

p0 <- ggplot(mf2[rt,], aes(x=Reliability,color=mqtlprobe)) + 
  geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p0,file="zhouvssugden_mqtl_retain.pdf",height=6,width=6)


p <- ggplot(mf2, aes(x=Reliability,color=MASK_general)) + 
  geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p,file="zhouvssugden.pdf",height=6,width=6)


p1 <- ggplot(mf2, aes(x=Reliability,color=MASK_sub30_copy)) + 
   geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p1,file="zhouvssugden_MASK_sub30_copy.pdf",height=6,width=6)

p2 <- ggplot(mf2, aes(x=Reliability,color=MASK_mapping)) + 
   geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p2,file="zhouvssugden_MASK_mapping.pdf",height=6,width=6)

p3 <- ggplot(mf2, aes(x=Reliability,color=MASK_extBase)) + 
   geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p3,file="zhouvssugden_MASK_extBase.pdf",height=6,width=6)

p4 <- ggplot(mf2, aes(x=Reliability,color=MASK_typeINextBaseSwitch)) + 
   geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p4,file="zhouvssugden_MASK_typeINextBaseSwitch.pdf",height=6,width=6)

p5 <- ggplot(mf2, aes(x=Reliability,color=MASK_snp5_GMAF1p)) + 
  geom_density() +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
  )

ggsave(p5,file="zhouvssugden_MASK.snp5.GMAF1p.pdf",height=6,width=6)

pdf("zhouvssugden_MASKcategories.pdf",height=6,width=8)
library(grid)
# Move to a new page
grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(p, vp = define_region(row = 1, col = 1))   # Span over two columns
print(p1, vp = define_region(row = 1, col = 2))
print(p2, vp = define_region(row = 2, col = 1))
print(p3, vp = define_region(row = 2, col = 2))
print(p4, vp = define_region(row = 3, col = 1))
print(p5, vp = define_region(row = 3, col = 2))


dev.off()





multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#MASK.mapping - whether the probe is masked for mapping reason. Probes retained should have high quality (>=40 on 0-60 scale) consistent (with designed MAPINFO) mapping (for both in the case of type I) without INDELs.
#MASK.typeINextBaseSwitch - whether the probe has a SNP in the extension base that causes a color channel switch from the official annotation (described as color-channel-switching, or CCS SNP in the reference). These probes should be processed differently than designed (by summing up both color channels instead of just the annotated color channel).
#MASK.rmsk15 - whether the 15bp 3'-subsequence of the probe overlap with repeat masker, this MASK is NOT recommended.
#MASK.sub25.copy, MASK.sub30.copy, MASK.sub35.copy and MASK.sub40.copy - whether the 25bp, 30bp, 35bp and 40bp 3'-subsequence of the probe is non-unique.
#MASK.snp5.common - whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the common SNPs from dbSNP (global MAF can be under 1%).
#MASK.snp5.GMAF1p - whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs with global MAF >1%.
#MASK.extBase - probes masked for extension base inconsistent with specified color channel (type-I) or CpG (type-II) based on mapping.
#MASK.general - recommended general purpose masking merged from "MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p".

