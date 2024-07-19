#Delfi Dorussen
#Import the DMR files for a mutant, join, and remove extra columns to save as a
#bed file

setwd("Z:/Delfi/MET1_WGBS/DMR_calling/")


mutants <- c("A_single", "B_single", "D_single", "AB_double", "AD_double", "BD_double","Aabbdd")

for(i in 1:7){
  mutant <- mutants[i]
  
  WT_mutant_DMRs_Chr1A <- read.table(paste("WT_",mutant,"_DMRs_Chr1A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr1B <- read.table(paste("WT_",mutant,"_DMRs_Chr1B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr1D <- read.table(paste("WT_",mutant,"_DMRs_Chr1D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr2A <- read.table(paste("WT_",mutant,"_DMRs_Chr2A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr2B <- read.table(paste("WT_",mutant,"_DMRs_Chr2B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr2D <- read.table(paste("WT_",mutant,"_DMRs_Chr2D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr3A <- read.table(paste("WT_",mutant,"_DMRs_Chr3A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr3B <- read.table(paste("WT_",mutant,"_DMRs_Chr3B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr3D <- read.table(paste("WT_",mutant,"_DMRs_Chr3D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr4A <- read.table(paste("WT_",mutant,"_DMRs_Chr4A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr4B <- read.table(paste("WT_",mutant,"_DMRs_Chr4B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr4D <- read.table(paste("WT_",mutant,"_DMRs_Chr4D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr5A <- read.table(paste("WT_",mutant,"_DMRs_Chr5A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr5B <- read.table(paste("WT_",mutant,"_DMRs_Chr5B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr5D <- read.table(paste("WT_",mutant,"_DMRs_Chr5D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr6A <- read.table(paste("WT_",mutant,"_DMRs_Chr6A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr6B <- read.table(paste("WT_",mutant,"_DMRs_Chr6B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr6D <- read.table(paste("WT_",mutant,"_DMRs_Chr6D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr7A <- read.table(paste("WT_",mutant,"_DMRs_Chr7A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr7B <- read.table(paste("WT_",mutant,"_DMRs_Chr7B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr7D <- read.table(paste("WT_",mutant,"_DMRs_Chr7D.txt",sep=""), header=TRUE)
  
  WT_mutant_DMRs <- rbind(WT_mutant_DMRs_Chr1A,WT_mutant_DMRs_Chr1B,WT_mutant_DMRs_Chr1D,
                          WT_mutant_DMRs_Chr2A,WT_mutant_DMRs_Chr2B,WT_mutant_DMRs_Chr2D,
                          WT_mutant_DMRs_Chr3A,WT_mutant_DMRs_Chr3B,WT_mutant_DMRs_Chr3D,
                          WT_mutant_DMRs_Chr4A,WT_mutant_DMRs_Chr4B,WT_mutant_DMRs_Chr4D,
                          WT_mutant_DMRs_Chr5A,WT_mutant_DMRs_Chr5B,WT_mutant_DMRs_Chr5D,
                          WT_mutant_DMRs_Chr6A,WT_mutant_DMRs_Chr6B,WT_mutant_DMRs_Chr6D,
                          WT_mutant_DMRs_Chr7A,WT_mutant_DMRs_Chr7B,WT_mutant_DMRs_Chr7D)
  
  head(WT_mutant_DMRs)
  
  WT_mutant_hypoDMRs <- subset(WT_mutant_DMRs, regionType=="loss")
  WT_mutant_hyperDMRs <- subset(WT_mutant_DMRs, regionType=="gain")
  
  WT_mutant_hypoDMRs_bed <- WT_mutant_hypoDMRs[,1:3]
  WT_mutant_hyperDMRs_bed <- WT_mutant_hyperDMRs[,1:3]
  
  write.table(WT_mutant_hypoDMRs_bed,paste(mutant,"_hypo_DMRs.bed",sep=""),col.names=FALSE, 
              row.names=FALSE, quote=FALSE, sep="\t")
  write.table(WT_mutant_hyperDMRs_bed,paste(mutant,"_hyper_DMRs.bed",sep=""),col.names=FALSE, 
              row.names=FALSE, quote=FALSE, sep="\t")
  
}

