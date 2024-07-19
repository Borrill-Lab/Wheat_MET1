#Delfi Dorussen
#Aim to call DMRs for each of the met1 mutants relative to the WT segregant
#Output: list of DMRs in the mutant as txt file

library(DMRcaller)

setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/")

chromosomes <- c("Chr1A","Chr1B","Chr1D","Chr2A","Chr2B","Chr2D","Chr3A","Chr3B","Chr3D",
                 "Chr4A","Chr4B","Chr4D","Chr5A","Chr5B","Chr5D","Chr6A","Chr6B","Chr6D",
                 "Chr7A","Chr7B","Chr7D")

for(i in 1:21){
  setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/")
  
  chromosome <- chromosomes[i]
  
  WT <- readBismark(paste("Unknown_BU703-001Z0001_1.",chromosome,".CX_report.txt", sep=""))
  mutant <- readBismark(paste("Unknown_BU703-001Z0007_1.",chromosome,".CX_report.txt", sep=""))
  
  methylationDataList <- list("WT"=WT, "mutant"=mutant)
  
  DMRsBinsCG <- computeDMRs(methylationDataList[["WT"]], methylationDataList[["mutant"]],
                            context="CG", method="bins", binSize=100, 
                            pValueThreshold=0.01, minCytosinesCount=10, minProportionDifference=0.4,
                            minGap=0, minReadsPerCytosine=4, cores=1)
  
  WT_mutant_DMRs <- as.data.frame(DMRsBinsCG)
  
  setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_calling/")
  
  write.table(WT_mutant_DMRs, file=paste("WT_BD_double_DMRs_",chromosome,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE)
}
