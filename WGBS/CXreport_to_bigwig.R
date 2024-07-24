#Delfi Dorussen
#Aim to split CX reports for each chromosome into the separate contexts (CG, CHG, CHH)
#and export each as a bigwig file, ultimately to use for plotting average methylation
#across genes and TEs

library(DMRcaller)
library(rtracklayer)

chromosomes <- c("Chr1A","Chr1B","Chr1D","Chr2A","Chr2B","Chr2D","Chr3A","Chr3B","Chr3D",
                 "Chr4A","Chr4B","Chr4D","Chr5A","Chr5B","Chr5D","Chr6A","Chr6B","Chr6D",
                 "Chr7A","Chr7B","Chr7D")

lengths <- c(594102053, 689851867, 495453184, 780798554, 801256712, 651852608, 750843637, 830829763, 615552420,
             744588156, 673617493, 509857065, 709773740, 713149756, 566080676, 618079258, 720988477, 473592717,
             736706231, 750620384, 638686045)

for(i in 7:9){
  
  for(j in 1:21){
    chromosome <- chromosomes[j]
    length <- lengths[j]
    
    setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/")
    
    CX_report <- readBismark(paste("Unknown_BU703-001Z000",i,"_1.",chromosome,".CX_report.txt",sep=""))
    
    CG_report <- CX_report[which(CX_report$context=="CG")]
    CHG_report <- CX_report[which(CX_report$context=="CHG")]
    CHH_report <- CX_report[which(CX_report$context=="CHH")]
    
    seqlengths(CG_report) <- length
    seqlengths(CHG_report) <- length
    seqlengths(CHH_report) <- length
    
    CG_report <- CG_report[CG_report$readsN >= 4]
    CHG_report <- CHG_report[CHG_report$readsN >= 4]
    CHH_report <- CHH_report[CHH_report$readsN >= 4]
    
    CG_report$score <- CG_report$readsM / CG_report$readsN
    CHG_report$score <- CHG_report$readsM / CHG_report$readsN
    CHH_report$score <- CHH_report$readsM / CHH_report$readsN
    
    setwd("/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/")
    
    rtracklayer::export.bw(CG_report, paste("Z000",i,"_",chromosome,"_CG.bw", sep=""))
    rtracklayer::export.bw(CHG_report, paste("Z000",i,"_",chromosome,"_CHG.bw", sep=""))
    rtracklayer::export.bw(CHH_report, paste("Z000",i,"_",chromosome,"_CHH.bw", sep=""))
  }
  
}

